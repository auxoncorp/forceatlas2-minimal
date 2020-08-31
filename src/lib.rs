#![feature(maybe_uninit_ref)]
#![feature(trait_alias)]

mod util;

pub use util::{clone_slice_mut, Coord, Edge, Node, PointIter, PointList, Position};

use itertools::izip;
#[cfg(feature = "rand")]
use rand::Rng;

#[derive(Clone)]
pub struct Settings<T: Coord> {
	pub dimensions: usize,
	pub kg: T,
	pub kr: T,
	pub lin_log: bool,
	pub prevent_overlapping: Option<T>,
	pub strong_gravity: bool,
}

impl<T: Coord> Default for Settings<T> {
	fn default() -> Self {
		Self {
			dimensions: 2,
			kg: T::one(),
			kr: T::one(),
			lin_log: false,
			prevent_overlapping: None,
			strong_gravity: false,
		}
	}
}

pub struct Layout<T: Coord> {
	pub edges: Vec<Edge>,
	pub nodes: Vec<Node>,
	pub points: PointList<T>,
	pub settings: Settings<T>,
	speeds: PointList<T>,
}

impl<'a, T: Coord> Layout<T> {
	#[cfg(feature = "rand")]
	pub fn from_graph(
		edges: Vec<Edge>,
		nb_nodes: usize,
		settings: Settings<T>,
		rand_range: (T, T),
	) -> Self
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
					(0..nb)
						.map(|_| rng.gen_range(rand_range.0.clone(), rand_range.1.clone()))
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
		if self.settings.lin_log {
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

				if let Some(node_size) = &self.settings.prevent_overlapping {
					let dprime = d.clone() - node_size.clone();
					if dprime.positive() {
						d = dprime;
					}
				}
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
		} else {
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

	fn apply_gravity(&mut self) {
		if self.settings.strong_gravity {
			for (n, pos, speed) in izip!(
				self.nodes.iter(),
				self.points.iter(),
				self.speeds.iter_mut()
			) {
				let f = util::gravity(self.settings.kg.clone(), n);
				for i in 0usize..self.settings.dimensions {
					speed[i] += pos[i].clone() - f.clone() * pos[i].clone();
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
				let f = util::gravity(self.settings.kg.clone(), n) / d;
				for i in 0usize..self.settings.dimensions {
					speed[i] += pos[i].clone() - f.clone() * pos[i].clone();
				}
			}
		}
	}

	fn apply_repulsion(&mut self) {
		for (n1, n2) in self.edges.iter() {
			let n1_pos = self.points.get(*n1);
			let mut d2 = T::zero();
			let mut di_v = self.points.get_clone(*n2);
			let di = di_v.as_mut_slice();
			for i in 0usize..self.settings.dimensions {
				di[i] -= n1_pos[i].clone();
				d2 += di[i].clone().pow_n(2u32);
			}
			if d2.is_zero() {
				continue;
			}

			let f = T::from(
				(unsafe { self.nodes.get_unchecked(*n1) }.degree + 1)
					* (unsafe { self.nodes.get_unchecked(*n2) }.degree + 1),
			) / d2 * self.settings.kg.clone();

			let n1_speed = self.speeds.get_mut(*n1);
			for i in 0usize..self.settings.dimensions {
				n1_speed[i] -= f.clone() * di[i].clone();
			}
			let n2_speed = self.speeds.get_mut(*n2);
			for i in 0usize..self.settings.dimensions {
				n2_speed[i] += f.clone() * di[i].clone();
			}
		}
	}

	fn apply_forces(&mut self) {
		for (pos, speed) in self.points.points.iter_mut().zip(self.speeds.points.iter()) {
			*pos += speed.clone();
		}
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_global() {
		let mut layout = Layout::<f64>::from_graph(
			vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 4)],
			5,
			Settings::default(),
			(-100.0, 100.0),
		);

		for _ in 0..10 {
			layout.iteration();
		}

		layout.points.iter().for_each(|pos| println!("{:?}", pos));
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
}
