#![feature(specialization)]
#![feature(trait_alias)]
#![allow(incomplete_features)]

mod util;

pub use util::{Coord, Edge, Node, PointIter, PointIterMut, PointList, Position};

use itertools::izip;

#[derive(Clone)]
pub struct Settings<T: Coord> {
	/// Number of spatial dimensions
	pub dimensions: usize,
	/// Move hubs (high degree nodes) to the center
	pub dissuade_hubs: bool,
	/// Attraction coefficient
	pub ka: T,
	/// Gravity coefficient
	pub kg: T,
	/// Repulsion coefficient
	pub kr: T,
	/// Logarithmic attraction
	pub lin_log: bool,
	/// Prevent node overlapping for a prettier graph (node_size, kr_prime).
	///
	/// `node_size` is the radius around a node where the repulsion coefficient is `kr_prime`.
	/// `kr_prime` is arbitrarily set to `100.0` in Gephi implementation.
	pub prevent_overlapping: Option<(T, T)>,
	/// Gravity does not decrease with distance, resulting in a more compact graph.
	pub strong_gravity: bool,
	/// Optimize repulsion using Barnes-Hut algorithm (time passes from N^2 to NlogN).
	/// The argument is theta.
	///
	/// **Note**: only implemented for `T=f64` and `dimension` 2 or 3.
	#[cfg(feature = "barnes_hut")]
	pub barnes_hut: Option<T>,
}

impl<T: Coord> Default for Settings<T> {
	fn default() -> Self {
		Self {
			dimensions: 2,
			dissuade_hubs: false,
			ka: T::from(0.5),
			kg: T::one(),
			kr: T::one(),
			lin_log: false,
			prevent_overlapping: None,
			strong_gravity: false,
			#[cfg(feature = "barnes_hut")]
			barnes_hut: None,
		}
	}
}

pub struct Layout<T: Coord> {
	pub edges: Vec<Edge>,
	pub nodes: Vec<Node>,
	/// List of the nodes' positions
	pub points: PointList<T>,
	pub settings: Settings<T>,
	speeds: PointList<T>,
}

impl<'a, T: Coord + std::fmt::Debug> Layout<T>
where
	Layout<T>: Repulsion,
{
	/// Instanciates a randomly positioned layout from a directed graph
	#[cfg(feature = "rand")]
	pub fn from_graph(edges: Vec<Edge>, nb_nodes: usize, settings: Settings<T>) -> Self
	where
		rand::distributions::Standard: rand::distributions::Distribution<T>,
		T: rand::distributions::uniform::SampleUniform,
	{
		let mut nodes: Vec<Node> = (0..nb_nodes).map(|_| Node { degree: 0 }).collect();
		for (n1, n2) in edges.iter() {
			nodes.get_mut(*n1).unwrap().degree += 1;
			nodes.get_mut(*n2).unwrap().degree += 1;
		}
		let nb = nb_nodes * settings.dimensions;
		Self {
			nodes,
			edges,
			points: PointList {
				dimensions: settings.dimensions,
				points: {
					let mut rng = rand::thread_rng();
					(0..nb_nodes)
						.flat_map(|_| util::sample_unit_nsphere(&mut rng, settings.dimensions))
						.collect()
				},
			},
			speeds: PointList {
				dimensions: settings.dimensions,
				points: (0..nb).map(|_| T::zero()).collect(),
			},
			settings,
		}
	}

	/// Instanciates layout from a directed graph, using initial positions
	pub fn from_position_graph<I: Iterator<Item = &'a Position<T>>>(
		edges: Vec<Edge>,
		nodes: I,
		settings: Settings<T>,
	) -> Self
	where
		T: 'a,
	{
		let mut points = Vec::new();
		let mut nodes: Vec<Node> = nodes
			.map(|pos| {
				points.extend_from_slice(pos);
				Node { degree: 0 }
			})
			.collect();

		for (n1, n2) in edges.iter() {
			nodes.get_mut(*n1).unwrap().degree += 1;
			nodes.get_mut(*n2).unwrap().degree += 1;
		}
		let nb = nodes.len() * settings.dimensions;
		Self {
			nodes,
			edges,
			points: PointList {
				dimensions: settings.dimensions,
				points,
			},
			speeds: PointList {
				dimensions: settings.dimensions,
				points: (0..nb).map(|_| T::zero()).collect(),
			},
			settings,
		}
	}

	/// Computes an iteration of ForceAtlas2
	pub fn iteration(&mut self) {
		self.init_iteration();
		self.apply_attraction();
		self.apply_repulsion();
		self.apply_gravity();
		self.apply_forces();
	}

	fn init_iteration(&mut self) {
		for speed in self.speeds.points.iter_mut() {
			*speed = T::zero();
		}
	}

	fn apply_attraction(&mut self) {
		if let Some((node_size, _)) = &self.settings.prevent_overlapping {
			if self.settings.lin_log {
				if self.settings.dissuade_hubs {
					// PO LOG DH
					for (n1, n2) in self.edges.iter() {
						let mut d = T::zero();
						let n1_pos = self.points.get(*n1);
						let mut di_v = self.points.get_clone(*n2);
						let di = di_v.as_mut_slice();
						for i in 0usize..self.settings.dimensions {
							di[i] -= n1_pos[i].clone();
							d += di[i].clone().pow_n(2u32);
						}
						d = d.sqrt();

						let dprime = d - node_size.clone();
						if dprime.non_positive() {
							continue;
						}
						let n1_degree = T::from(self.nodes.get(*n1).unwrap().degree);
						let f =
							dprime.clone().ln_1p() / dprime / n1_degree * self.settings.ka.clone();

						let n1_speed = self.speeds.get_mut(*n1);
						for i in 0usize..self.settings.dimensions {
							n1_speed[i] += f.clone() * di[i].clone();
						}
						let n2_speed = self.speeds.get_mut(*n2);
						for i in 0usize..self.settings.dimensions {
							n2_speed[i] -= f.clone() * di[i].clone();
						}
					}
				} else {
					// PO LOG --
					for (n1, n2) in self.edges.iter() {
						let mut d = T::zero();
						let n1_pos = self.points.get(*n1);
						let mut di_v = self.points.get_clone(*n2);
						let di = di_v.as_mut_slice();
						for i in 0usize..self.settings.dimensions {
							di[i] -= n1_pos[i].clone();
							d += di[i].clone().pow_n(2u32);
						}
						d = d.sqrt();

						let dprime = d - node_size.clone();
						if dprime.non_positive() {
							continue;
						}
						let f = dprime.clone().ln_1p() / dprime;

						let n1_speed = self.speeds.get_mut(*n1);
						for i in 0usize..self.settings.dimensions {
							n1_speed[i] += f.clone() * di[i].clone();
						}
						let n2_speed = self.speeds.get_mut(*n2);
						for i in 0usize..self.settings.dimensions {
							n2_speed[i] -= f.clone() * di[i].clone();
						}
					}
				}
			} else {
				if self.settings.dissuade_hubs {
					// PO --- DH
					for (n1, n2) in self.edges.iter() {
						let mut d = T::zero();
						let n1_pos = self.points.get(*n1);
						let mut di_v = self.points.get_clone(*n2);
						let di = di_v.as_mut_slice();
						for i in 0usize..self.settings.dimensions {
							di[i] -= n1_pos[i].clone();
							d += di[i].clone().pow_n(2u32);
						}
						d = d.sqrt();

						let dprime = d.clone() - node_size.clone();
						if dprime.non_positive() {
							dbg!(dprime);
							continue;
						}
						let n1_degree = T::from(self.nodes.get(*n1).unwrap().degree);
						let f = dprime / d / n1_degree;

						let n1_speed = self.speeds.get_mut(*n1);
						for i in 0usize..self.settings.dimensions {
							n1_speed[i] += f.clone() * di[i].clone();
						}
						let n2_speed = self.speeds.get_mut(*n2);
						for i in 0usize..self.settings.dimensions {
							n2_speed[i] -= f.clone() * di[i].clone();
						}
					}
				} else {
					// PO --- --
					for (n1, n2) in self.edges.iter() {
						let mut d = T::zero();
						let n1_pos = self.points.get(*n1).clone();
						let mut di_v = self.points.get_clone(*n2);
						let di = di_v.as_mut_slice();
						for i in 0usize..self.settings.dimensions {
							di[i] -= n1_pos[i].clone();
							d += di[i].clone().pow_n(2u32);
						}
						d = d.sqrt();

						let dprime = d.clone() - node_size.clone();
						if dprime.non_positive() {
							continue;
						}
						let f = dprime / d;

						let n1_speed = self.speeds.get_mut(*n1);
						for i in 0usize..self.settings.dimensions {
							n1_speed[i] += f.clone() * di[i].clone();
						}
						let n2_speed = self.speeds.get_mut(*n2);
						for i in 0usize..self.settings.dimensions {
							n2_speed[i] -= f.clone() * di[i].clone();
						}
					}
				}
			}
		} else {
			if self.settings.lin_log {
				if self.settings.dissuade_hubs {
					// -- LOG DH
					for (n1, n2) in self.edges.iter() {
						let mut d = T::zero();
						let mut di_v = self.points.get_clone(*n2);
						let di = di_v.as_mut_slice();
						for i in 0usize..self.settings.dimensions {
							di[i] -= self.points.get(*n1)[i].clone();
							d += di[i].clone().pow_n(2u32);
						}
						if d.is_zero() {
							continue;
						}
						d = d.sqrt();

						let n1_degree = T::from(self.nodes.get(*n1).unwrap().degree);
						let f = d.clone().ln_1p() / d / n1_degree;

						let n1_speed = self.speeds.get_mut(*n1);
						for i in 0usize..self.settings.dimensions {
							n1_speed[i] += f.clone() * di[i].clone();
						}
						let n2_speed = self.speeds.get_mut(*n2);
						for i in 0usize..self.settings.dimensions {
							n2_speed[i] -= f.clone() * di[i].clone();
						}
					}
				} else {
					// -- LOG --
					for (n1, n2) in self.edges.iter() {
						let mut d = T::zero();
						let mut di_v = self.points.get_clone(*n2);
						let di = di_v.as_mut_slice();
						for i in 0usize..self.settings.dimensions {
							di[i] -= self.points.get(*n1)[i].clone();
							d += di[i].clone().pow_n(2u32);
						}
						if d.is_zero() {
							continue;
						}
						d = d.sqrt();

						let f = d.clone().ln_1p() / d;

						let n1_speed = self.speeds.get_mut(*n1);
						for i in 0usize..self.settings.dimensions {
							n1_speed[i] += f.clone() * di[i].clone();
						}
						let n2_speed = self.speeds.get_mut(*n2);
						for i in 0usize..self.settings.dimensions {
							n2_speed[i] -= f.clone() * di[i].clone();
						}
					}
				}
			} else {
				if self.settings.dissuade_hubs {
					// -- --- DH
					for (n1, n2) in self.edges.iter() {
						let n1_speed = self.speeds.get_mut(*n1);
						let n1_pos = self.points.get(*n1).clone();
						let mut di_v = self.points.get_clone(*n2);
						let di = di_v.as_mut_slice();
						let n1_degree = T::from(self.nodes.get(*n1).unwrap().degree);
						for i in 0usize..self.settings.dimensions {
							di[i] -= n1_pos[i].clone();
							di[i] /= n1_degree.clone();
							n1_speed[i] += di[i].clone();
						}
						let n2_speed = self.speeds.get_mut(*n2);
						for i in 0usize..self.settings.dimensions {
							n2_speed[i] -= di[i].clone();
						}
					}
				} else {
					// -- --- --
					for (n1, n2) in self.edges.iter() {
						let n1_speed = self.speeds.get_mut(*n1);
						let n1_pos = self.points.get(*n1).clone();
						let mut di_v = self.points.get_clone(*n2);
						let di = di_v.as_mut_slice();
						for i in 0usize..self.settings.dimensions {
							di[i] -= n1_pos[i].clone();
							n1_speed[i] += di[i].clone();
						}
						let n2_speed = self.speeds.get_mut(*n2);
						for i in 0usize..self.settings.dimensions {
							n2_speed[i] -= di[i].clone();
						}
					}
				}
			}
		}
	}

	fn apply_gravity(&mut self) {
		if self.settings.strong_gravity {
			for (n, pos, speed) in izip!(
				self.nodes.iter(),
				self.points.iter(),
				self.speeds.iter_mut()
			) {
				let f = T::from(n.degree + 1) * self.settings.kg.clone();
				for i in 0usize..self.settings.dimensions {
					speed[i] -= f.clone() * pos[i].clone();
				}
			}
		} else {
			for (n, pos, speed) in izip!(
				self.nodes.iter(),
				self.points.iter(),
				self.speeds.iter_mut()
			) {
				let d = util::norm(pos);
				if d.is_zero() {
					continue;
				}
				let f = T::from(n.degree + 1) * self.settings.kg.clone() / d;
				for i in 0usize..self.settings.dimensions {
					speed[i] -= f.clone() * pos[i].clone();
				}
			}
		}
	}

	#[inline]
	fn inner_apply_repulsion(&mut self) {
		if let Some((node_size, krprime)) = &self.settings.prevent_overlapping {
			for n1 in 0..self.nodes.len() - 1 {
				let n1_pos = self.points.get(n1);
				for n2 in n1 + 1..self.nodes.len() {
					let mut d = T::zero();
					let mut di_v = self.points.get_clone(n2);
					let di = di_v.as_mut_slice();
					for i in 0usize..self.settings.dimensions {
						di[i] -= n1_pos[i].clone();
						d += di[i].clone().pow_n(2u32);
					}
					d = d.sqrt();

					let dprime = d.clone() - node_size.clone();
					let f = if dprime.positive() {
						self.settings.kr.clone() / dprime
					} else if dprime.is_zero() {
						continue;
					} else {
						krprime.clone()
					} * T::from(
						(unsafe { self.nodes.get_unchecked(n1) }.degree + 1)
							* (unsafe { self.nodes.get_unchecked(n2) }.degree + 1),
					) / d;

					let n1_speed = self.speeds.get_mut(n1);
					for i in 0usize..self.settings.dimensions {
						n1_speed[i] -= f.clone() * di[i].clone();
					}
					let n2_speed = self.speeds.get_mut(n2);
					for i in 0usize..self.settings.dimensions {
						n2_speed[i] += f.clone() * di[i].clone();
					}
				}
			}
		} else {
			for n1 in 0..self.nodes.len() - 1 {
				let n1_pos = self.points.get(n1);
				for n2 in n1 + 1..self.nodes.len() {
					let mut d2 = T::zero();
					let mut di_v = self.points.get_clone(n2);
					let di = di_v.as_mut_slice();
					for i in 0usize..self.settings.dimensions {
						di[i] -= n1_pos[i].clone();
						d2 += di[i].clone().pow_n(2u32);
					}
					if d2.is_zero() {
						continue;
					}
					let d = d2.sqrt();

					let f = T::from(
						(unsafe { self.nodes.get_unchecked(n1) }.degree + 1)
							* (unsafe { self.nodes.get_unchecked(n2) }.degree + 1),
					) / d * self.settings.kr.clone();

					let n1_speed = self.speeds.get_mut(n1);
					for i in 0usize..self.settings.dimensions {
						n1_speed[i] -= f.clone() * di[i].clone();
					}
					let n2_speed = self.speeds.get_mut(n2);
					for i in 0usize..self.settings.dimensions {
						n2_speed[i] += f.clone() * di[i].clone();
					}
				}
			}
		}
	}

	fn apply_forces(&mut self) {
		for (pos, speed) in self.points.points.iter_mut().zip(self.speeds.points.iter()) {
			*pos += speed.clone();
		}
	}
}

#[cfg(feature = "barnes_hut")]
impl Layout<f64> {
	#[inline]
	fn inner_apply_repulsion_barnes_hut(&mut self) {
		match self.settings.dimensions {
			2 => {
				let particles: Vec<nbody_barnes_hut::particle_2d::Particle2D> = self
					.points
					.iter()
					.zip(self.nodes.iter())
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
					self.settings.barnes_hut.unwrap(),
				);
				let kr = self.settings.kr.clone();
				if let Some((node_size, krprime)) = &self.settings.prevent_overlapping {
					izip!(
						particles.into_iter(),
						self.speeds.iter_mut(),
						self.nodes.iter()
					)
					.for_each(|(particle, speed, node)| {
						let nbody_barnes_hut::vector_2d::Vector2D { x, y } = tree
							.calc_forces_on_particle(
								particle.position,
								node.degree + 1,
								|d2, m1, dv, m2| {
									let d = d2.sqrt();
									let dprime = d - node_size;
									(if dprime.positive() {
										kr / dprime
									} else if dprime.is_zero() {
										return nbody_barnes_hut::vector_2d::Vector2D {
											x: 0.0,
											y: 0.0,
										};
									} else {
										*krprime
									}) * m1 * m2 as f64 / d * dv
								},
							);
						speed[0] -= x;
						speed[1] -= y;
					});
				} else {
					izip!(
						particles.into_iter(),
						self.speeds.iter_mut(),
						self.nodes.iter()
					)
					.for_each(|(particle, speed, node)| {
						let nbody_barnes_hut::vector_2d::Vector2D { x, y } = tree
							.calc_forces_on_particle(
								particle.position,
								node.degree + 1,
								|d2, m1, dv, m2| m2 as f64 * m1 / d2.sqrt() * kr * dv,
							);
						speed[0] -= x;
						speed[1] -= y;
					});
				}
			}
			3 => {
				let particles: Vec<nbody_barnes_hut::particle_3d::Particle3D> = self
					.points
					.iter()
					.zip(self.nodes.iter())
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
					self.settings.barnes_hut.unwrap(),
				);
				let kr = self.settings.kr.clone();
				if let Some((node_size, krprime)) = &self.settings.prevent_overlapping {
					izip!(
						particles.into_iter(),
						self.speeds.iter_mut(),
						self.nodes.iter()
					)
					.for_each(|(particle, speed, node)| {
						let nbody_barnes_hut::vector_3d::Vector3D { x, y, z } = tree
							.calc_forces_on_particle(
								particle.position,
								node.degree + 1,
								|d2, m1, dv, m2| {
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
										*krprime
									}) * m1 * m2 as f64 / d * dv
								},
							);
						speed[0] -= x;
						speed[1] -= y;
						speed[2] -= z;
					});
				} else {
					izip!(
						particles.into_iter(),
						self.speeds.iter_mut(),
						self.nodes.iter()
					)
					.for_each(|(particle, speed, node)| {
						let nbody_barnes_hut::vector_3d::Vector3D { x, y, z } = tree
							.calc_forces_on_particle(
								particle.position,
								node.degree + 1,
								|d2, m1, dv, m2| m2 as f64 * m1 / d2.sqrt() * kr * dv,
							);
						speed[0] -= x;
						speed[1] -= y;
						speed[2] -= z;
					});
				}
			}
			_ => unimplemented!("Barnes-Hut only implemented for 2D and 3D"),
		}
	}
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
}

#[cfg(test)]
mod tests {
	use super::*;

	#[cfg(feature = "rand")]
	#[test]
	fn test_global() {
		let mut layout = Layout::<f64>::from_graph(
			vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 4)],
			5,
			Settings::default(),
		);

		for _ in 0..10 {
			layout.iteration();
		}

		layout.points.iter().for_each(|pos| println!("{:?}", pos));
	}

	#[test]
	fn test_apply_forces() {
		let mut layout = Layout::<f64>::from_position_graph(
			vec![(0, 1)],
			vec![vec![-1.0, -1.0].as_slice(), vec![1.0, 1.0].as_slice()].into_iter(),
			Settings::default(),
		);
		layout
			.speeds
			.points
			.iter_mut()
			.enumerate()
			.for_each(|(i, s)| *s += i as f64);
		layout.apply_forces();
		assert_eq!(layout.points.points, vec![-1.0, 0.0, 3.0, 4.0]);
	}

	#[test]
	fn test_init_iteration() {
		let mut layout = Layout::<f64>::from_position_graph(
			vec![(0, 1)],
			vec![vec![-1.0, -1.0].as_slice(), vec![1.0, 1.0].as_slice()].into_iter(),
			Settings::default(),
		);
		layout
			.speeds
			.points
			.iter_mut()
			.enumerate()
			.for_each(|(i, s)| *s += i as f64);
		layout.init_iteration();
		assert_eq!(layout.speeds.points, vec![0.0, 0.0, 0.0, 0.0]);
	}

	#[test]
	fn test_forces() {
		let mut layout = Layout::<f64>::from_position_graph(
			vec![(0, 1)],
			vec![vec![-1.0, -1.0].as_slice(), vec![1.0, 1.0].as_slice()].into_iter(),
			Settings::default(),
		);

		layout.init_iteration();
		layout.apply_attraction();

		let speed_1 = dbg!(layout.speeds.get(0));
		let speed_2 = dbg!(layout.speeds.get(1));

		assert!(speed_1[0] > 0.0);
		assert!(speed_1[1] > 0.0);
		assert!(speed_2[0] < 0.0);
		assert!(speed_2[1] < 0.0);

		layout.init_iteration();
		layout.apply_repulsion();

		let speed_1 = dbg!(layout.speeds.get(0));
		let speed_2 = dbg!(layout.speeds.get(1));

		assert!(speed_1[0] < 0.0);
		assert!(speed_1[1] < 0.0);
		assert!(speed_2[0] > 0.0);
		assert!(speed_2[1] > 0.0);

		layout.init_iteration();
		layout.apply_gravity();

		let speed_1 = dbg!(layout.speeds.get(0));
		let speed_2 = dbg!(layout.speeds.get(1));

		assert!(speed_1[0] > 0.0);
		assert!(speed_1[1] > 0.0);
		assert!(speed_2[0] < 0.0);
		assert!(speed_2[1] < 0.0);
	}

	#[cfg(feature = "barnes_hut")]
	#[test]
	fn test_barnes_hut_2d() {
		let mut layout = Layout::<f64>::from_position_graph(
			vec![(0, 1)],
			vec![vec![-1.0, -1.0].as_slice(), vec![1.0, 1.0].as_slice()].into_iter(),
			Settings::default(),
		);

		layout.settings.barnes_hut = Some(0.5);
		layout.init_iteration();
		layout.apply_repulsion();

		let speed_1 = dbg!(layout.speeds.get(0));
		let speed_2 = dbg!(layout.speeds.get(1));

		assert!(speed_1[0] < 0.0);
		assert!(speed_1[1] < 0.0);
		assert!(speed_2[0] > 0.0);
		assert!(speed_2[1] > 0.0);
	}

	#[test]
	fn test_convergence() {
		let mut layout = Layout::<f64>::from_position_graph(
			vec![(0, 1), (1, 2)],
			vec![
				vec![-1.1, -1.0].as_slice(),
				vec![0.0, 0.0].as_slice(),
				vec![1.0, 1.0].as_slice(),
			]
			.into_iter(),
			Settings {
				dimensions: 2,
				dissuade_hubs: false,
				ka: 0.5,
				kg: 0.01,
				kr: 0.01,
				lin_log: false,
				prevent_overlapping: None,
				strong_gravity: false,
				#[cfg(feature = "barnes_hut")]
				barnes_hut: None,
			},
		);

		for _ in 0..10 {
			println!("new iteration");
			layout.init_iteration();
			layout.apply_attraction();
			println!("{:?}", layout.speeds.points);
			layout.init_iteration();
			layout.apply_repulsion();
			println!("{:?}", layout.speeds.points);
			layout.init_iteration();
			layout.apply_gravity();
			println!("{:?}", layout.speeds.points);
			layout.apply_forces();
			//layout.iteration();

			dbg!(&layout.points.points);
			let point_1 = layout.points.get(0);
			let point_2 = layout.points.get(1);
			dbg!(((point_2[0] - point_1[0]).powi(2) + (point_2[1] - point_1[1]).powi(2)).sqrt());
		}
	}
}
