use crate::{layout::Layout, util::*};

use itertools::izip;

pub fn apply_repulsion<T: Coord + std::fmt::Debug>(layout: &mut Layout<T>) {
	let mut di = valloc(layout.settings.dimensions);
	let mut n2_iter = layout.points.iter();
	for (n1, (n1_node, n1_pos)) in layout.nodes.iter().zip(layout.points.iter()).enumerate() {
		let n1_degree = n1_node.degree + 1;
		n2_iter.offset = (n1 + 1) * layout.settings.dimensions;
		for (n2, n2_pos) in (n1 + 1..).zip(&mut n2_iter) {
			di.clone_from_slice(n2_pos);

			let d2 = di
				.iter_mut()
				.zip(n1_pos.iter())
				.map(|(di, n1_pos)| {
					*di -= n1_pos.clone();
					di.clone().pow_n(2u32)
				})
				.sum::<T>();
			if d2.is_zero() {
				continue;
			}

			let f = T::from(n1_degree * (unsafe { layout.nodes.get_unchecked(n2) }.degree + 1))
				/ d2 * layout.settings.kr.clone();

			let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);
			izip!(n1_speed.iter_mut(), n2_speed.iter_mut(), di.iter()).for_each(
				|(n1_speed, n2_speed, di)| {
					let s = f.clone() * di.clone();
					*n1_speed -= s.clone();
					*n2_speed += s;
				},
			);
		}
	}
}

pub fn apply_repulsion_fast<T: Copy + Coord + std::fmt::Debug>(layout: &mut Layout<T>) {
	let mut di = valloc(layout.settings.dimensions);
	let mut n2_iter = layout.points.iter();
	let last_dim = layout.settings.dimensions - 1;
	for (n1, (n1_node, n1_pos)) in layout.nodes.iter().zip(layout.points.iter()).enumerate() {
		let n1_degree = n1_node.degree + 1;
		n2_iter.offset = (n1 + 1) * layout.settings.dimensions;
		for (n2, n2_pos) in (n1 + 1..).zip(&mut n2_iter) {
			di.copy_from_slice(n2_pos);

			/*let mut i = 0usize;
			let mut d2 = T::zero();
			loop {
				let di = unsafe{di.get_unchecked_mut(i)};
				let n1_pos = *unsafe{n1_pos.get_unchecked(i)};
				*di = unsafe{std::intrinsics::fsub_fast(*di, n1_pos)};
				d2 = unsafe{std::intrinsics::fadd_fast(d2, std::intrinsics::fmul_fast(*di, *di))};

				if i == last_dim {
					break;
				}
				i += 1;
			}*/

			let d2 = di
				.iter_mut()
				.zip(n1_pos.iter())
				.fold(T::zero(), |sum, (di, n1_pos)| {
					*di = unsafe { std::intrinsics::fsub_fast(*di, *n1_pos) };
					unsafe { std::intrinsics::fadd_fast(sum, std::intrinsics::fmul_fast(*di, *di)) }
				});

			if d2.is_zero() {
				continue;
			}

			let f = unsafe {
				std::intrinsics::fmul_fast(
					std::intrinsics::fdiv_fast(
						T::from(n1_degree * (layout.nodes.get_unchecked(n2).degree + 1)),
						d2,
					),
					layout.settings.kr,
				)
			};

			let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);
			let mut i = 0usize;
			loop {
				let di = unsafe { di.get_unchecked(i) };
				let n1_speed = unsafe { n1_speed.get_unchecked_mut(i) };
				let n2_speed = unsafe { n2_speed.get_unchecked_mut(i) };

				let s = unsafe { std::intrinsics::fmul_fast(f, *di) };
				*n1_speed = unsafe { std::intrinsics::fadd_fast(*n1_speed, s) };
				*n2_speed = unsafe { std::intrinsics::fsub_fast(*n2_speed, s) };

				if i == last_dim {
					break;
				}
				i += 1;
			}

			/*izip!(n1_speed.iter_mut(), n2_speed.iter_mut(), di.iter()).for_each(
				|(n1_speed, n2_speed, di)| {
					let s = unsafe{std::intrinsics::fmul_fast(f, *di)};
					*n1_speed -= s;
					*n2_speed += s;
				},
			);*/
		}
	}
}

pub fn apply_repulsion_fast_2d<T: Copy + Coord + std::fmt::Debug>(layout: &mut Layout<T>) {
	let mut n2_iter = layout.points.iter();
	for (n1, (n1_node, n1_pos)) in layout.nodes.iter().zip(layout.points.iter()).enumerate() {
		let n1_degree = n1_node.degree + 1;
		for (n2, n2_pos) in (n1 + 1..).zip(&mut n2_iter) {
			let dx = unsafe { *n2_pos.get_unchecked(0) - *n1_pos.get_unchecked(0) };
			let dy = unsafe { *n2_pos.get_unchecked(1) - *n1_pos.get_unchecked(1) };

			let d2 = dx * dx + dy * dy;
			if d2.is_zero() {
				continue;
			}

			let f = T::from(n1_degree * (unsafe { layout.nodes.get_unchecked(n2) }.degree + 1))
				/ d2 * layout.settings.kr;

			let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);
			unsafe { n1_speed.get_unchecked_mut(0) }.add_assign(f * dx);
			unsafe { n1_speed.get_unchecked_mut(1) }.add_assign(f * dy);
			unsafe { n2_speed.get_unchecked_mut(0) }.sub_assign(f * dx);
			unsafe { n2_speed.get_unchecked_mut(1) }.sub_assign(f * dy);
		}
	}
}

pub fn apply_repulsion_po<T: Coord + std::fmt::Debug>(layout: &mut Layout<T>) {
	let mut di = valloc(layout.settings.dimensions);
	let mut n2_iter = layout.points.iter();
	let (node_size, krprime) = unsafe {
		layout
			.settings
			.prevent_overlapping
			.as_ref()
			.unwrap_unchecked()
	};
	for (n1, (n1_node, n1_pos)) in layout.nodes.iter().zip(layout.points.iter()).enumerate() {
		let n1_degree = n1_node.degree + 1;
		n2_iter.offset = (n1 + 1) * layout.settings.dimensions;
		for (n2, n2_pos) in (n1 + 1..).zip(&mut n2_iter) {
			di.clone_from_slice(n2_pos);

			let d2 = di
				.iter_mut()
				.zip(n1_pos.iter())
				.map(|(di, n1_pos)| {
					*di -= n1_pos.clone();
					di.clone().pow_n(2u32)
				})
				.sum::<T>();
			if d2.is_zero() {
				continue;
			}

			let d = d2.clone().sqrt();
			let dprime = d.clone() - node_size.clone();

			let f = T::from(n1_degree * (unsafe { layout.nodes.get_unchecked(n2) }.degree + 1))
				/ d2 * if dprime.positive() {
				layout.settings.kr.clone() / dprime
			} else {
				krprime.clone()
			};

			let (n1_speed, n2_speed) = layout.speeds.get_2_mut(n1, n2);
			izip!(n1_speed.iter_mut(), n2_speed.iter_mut(), di.iter()).for_each(
				|(n1_speed, n2_speed, di)| {
					let s = f.clone() * di.clone();
					*n1_speed -= s.clone();
					*n2_speed += s;
				},
			);
		}
	}
}

#[cfg(feature = "barnes_hut")]
pub fn apply_repulsion_bh_2d(layout: &mut Layout<f64>) {
	let particles: Vec<nbody_barnes_hut::particle_2d::Particle2D> = layout
		.points
		.iter()
		.zip(layout.nodes.iter())
		.map(|(point, node)| nbody_barnes_hut::particle_2d::Particle2D {
			position: nbody_barnes_hut::vector_2d::Vector2D {
				x: point[0],
				y: point[1],
			},
			mass: (node.degree + 1) as f64,
		})
		.collect();
	let tree = nbody_barnes_hut::barnes_hut_2d::QuadTree::new(
		&particles
			.iter()
			.collect::<Vec<&nbody_barnes_hut::particle_2d::Particle2D>>(),
		layout.settings.barnes_hut.unwrap(),
	);
	let kr = layout.settings.kr;

	izip!(
		particles.into_iter(),
		layout.speeds.iter_mut(),
		layout.nodes.iter()
	)
	.for_each(|(particle, speed, node)| {
		let nbody_barnes_hut::vector_2d::Vector2D { x, y } =
			tree.calc_forces_on_particle(particle.position, node.degree + 1, |d2, m1, dv, m2| {
				unsafe {
					std::intrinsics::fmul_fast(
						std::intrinsics::fmul_fast(
							m2 as f64,
							std::intrinsics::fdiv_fast(m1, d2.sqrt()),
						),
						kr,
					) * dv
				}
				//m2 as f64 * m1 / d2.sqrt() * kr * dv
			});
		speed[0] -= x;
		speed[1] -= y;
	});
}

#[cfg(feature = "barnes_hut")]
pub fn apply_repulsion_bh_2d_po(layout: &mut Layout<f64>) {
	let particles: Vec<nbody_barnes_hut::particle_2d::Particle2D> = layout
		.points
		.iter()
		.zip(layout.nodes.iter())
		.map(|(point, node)| nbody_barnes_hut::particle_2d::Particle2D {
			position: nbody_barnes_hut::vector_2d::Vector2D {
				x: point[0],
				y: point[1],
			},
			mass: (node.degree + 1) as f64,
		})
		.collect();
	let tree = nbody_barnes_hut::barnes_hut_2d::QuadTree::new(
		&particles
			.iter()
			.collect::<Vec<&nbody_barnes_hut::particle_2d::Particle2D>>(),
		layout.settings.barnes_hut.unwrap(),
	);
	let kr = layout.settings.kr;
	let (node_size, krprime) = unsafe { layout.settings.prevent_overlapping.unwrap_unchecked() };
	izip!(
		particles.into_iter(),
		layout.speeds.iter_mut(),
		layout.nodes.iter()
	)
	.for_each(|(particle, speed, node)| {
		let nbody_barnes_hut::vector_2d::Vector2D { x, y } =
			tree.calc_forces_on_particle(particle.position, node.degree + 1, |d2, m1, dv, m2| {
				let d = d2.sqrt();
				let dprime = d - node_size;
				(if dprime.positive() {
					kr / dprime
				} else if dprime.is_zero() {
					return nbody_barnes_hut::vector_2d::Vector2D { x: 0.0, y: 0.0 };
				} else {
					krprime
				}) * m1 * m2 as f64 / d
					* dv
			});
		speed[0] -= x;
		speed[1] -= y;
	});
}

#[cfg(feature = "barnes_hut")]
pub fn apply_repulsion_bh_3d(layout: &mut Layout<f64>) {
	let particles: Vec<nbody_barnes_hut::particle_3d::Particle3D> = layout
		.points
		.iter()
		.zip(layout.nodes.iter())
		.map(|(point, node)| nbody_barnes_hut::particle_3d::Particle3D {
			position: nbody_barnes_hut::vector_3d::Vector3D {
				x: point[0],
				y: point[1],
				z: point[2],
			},
			mass: (node.degree + 1) as f64,
		})
		.collect();
	let tree = nbody_barnes_hut::barnes_hut_3d::OctTree::new(
		&particles
			.iter()
			.collect::<Vec<&nbody_barnes_hut::particle_3d::Particle3D>>(),
		layout.settings.barnes_hut.unwrap(),
	);
	let kr = layout.settings.kr;
	izip!(
		particles.into_iter(),
		layout.speeds.iter_mut(),
		layout.nodes.iter()
	)
	.for_each(|(particle, speed, node)| {
		let nbody_barnes_hut::vector_3d::Vector3D { x, y, z } =
			tree.calc_forces_on_particle(particle.position, node.degree + 1, |d2, m1, dv, m2| {
				m2 as f64 * m1 / d2.sqrt() * kr * dv
			});
		speed[0] -= x;
		speed[1] -= y;
		speed[2] -= z;
	});
}

#[cfg(feature = "barnes_hut")]
pub fn apply_repulsion_bh_3d_po(layout: &mut Layout<f64>) {
	let particles: Vec<nbody_barnes_hut::particle_3d::Particle3D> = layout
		.points
		.iter()
		.zip(layout.nodes.iter())
		.map(|(point, node)| nbody_barnes_hut::particle_3d::Particle3D {
			position: nbody_barnes_hut::vector_3d::Vector3D {
				x: point[0],
				y: point[1],
				z: point[2],
			},
			mass: (node.degree + 1) as f64,
		})
		.collect();
	let tree = nbody_barnes_hut::barnes_hut_3d::OctTree::new(
		&particles
			.iter()
			.collect::<Vec<&nbody_barnes_hut::particle_3d::Particle3D>>(),
		layout.settings.barnes_hut.unwrap(),
	);
	let kr = layout.settings.kr;
	let (node_size, krprime) = unsafe { layout.settings.prevent_overlapping.unwrap_unchecked() };
	izip!(
		particles.into_iter(),
		layout.speeds.iter_mut(),
		layout.nodes.iter()
	)
	.for_each(|(particle, speed, node)| {
		let nbody_barnes_hut::vector_3d::Vector3D { x, y, z } =
			tree.calc_forces_on_particle(particle.position, node.degree + 1, |d2, m1, dv, m2| {
				let d = d2.sqrt();
				let dprime = d - node_size;
				(if dprime.positive() {
					kr / dprime
				} else if dprime.is_zero() {
					return nbody_barnes_hut::vector_3d::Vector3D {
						x: 0.0,
						y: 0.0,
						z: 0.0,
					};
				} else {
					krprime
				}) * m1 * m2 as f64 / d
					* dv
			});
		speed[0] -= x;
		speed[1] -= y;
		speed[2] -= z;
	});
}

/*#[cfg(feature = "barnes_hut")]
impl Layout<f64> {

}

// This trait was needed to allow specialization
#[doc(hidden)]
pub trait Repulsion {
	fn apply_repulsion(&mut self);
}

#[cfg(not(feature = "barnes_hut"))]
impl<T: Coord + std::fmt::Debug> Repulsion for Layout<T> {
	fn apply_repulsion(&mut self) {
		self.inner_apply_repulsion()
	}
}

#[cfg(feature = "barnes_hut")]
default impl<T: Coord + std::fmt::Debug> Repulsion for Layout<T> {
	fn apply_repulsion(&mut self) {
		if self.settings.barnes_hut.is_none() {
			self.inner_apply_repulsion()
		} else {
			unimplemented!("Barnes-Hut only implemented for Layout<f64>")
		}
	}
}

#[cfg(feature = "barnes_hut")]
impl Repulsion for Layout<f64> {
	fn apply_repulsion(&mut self) {
		if self.settings.barnes_hut.is_some() {
			self.inner_apply_repulsion_barnes_hut()
		} else {
			self.inner_apply_repulsion()
		}
	}
}*/