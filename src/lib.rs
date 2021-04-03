#![feature(specialization)]
#![feature(trait_alias)]
#![allow(incomplete_features)]
#![feature(option_result_unwrap_unchecked)]
#![feature(core_intrinsics)]

mod forces;
mod layout;
mod util;

use forces::Repulsion;

pub use layout::{Layout, Settings};
pub use util::{Coord, Edge, Node, PointIter, PointIterMut, PointList, Position};

use itertools::izip;

impl<'a, T: Coord + std::fmt::Debug> Layout<T>
where
	Layout<T>: forces::Repulsion<T>,
{
	/// Instanciates a randomly positioned layout from an undirected graph
	///
	/// Assumes edges `(n1, n2)` respect `n1 < n2`
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
						.flat_map(|_| util::sample_unit_ncube(&mut rng, settings.dimensions))
						.collect()
				},
			},
			speed: T::one(),
			speeds: PointList {
				dimensions: settings.dimensions,
				points: (0..nb).map(|_| T::zero()).collect(),
			},
			old_speeds: PointList {
				dimensions: settings.dimensions,
				points: (0..nb).map(|_| T::zero()).collect(),
			},
			fn_attraction: forces::choose_attraction(&settings),
			fn_gravity: forces::choose_gravity(&settings),
			fn_repulsion: Self::choose_repulsion(&settings),
			settings,
		}
	}

	/// Instanciates layout from an undirected graph, using initial positions
	///
	/// Assumes edges `(n1, n2)` respect `n1 < n2`
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
			speed: T::one(),
			speeds: PointList {
				dimensions: settings.dimensions,
				points: (0..nb).map(|_| T::zero()).collect(),
			},
			old_speeds: PointList {
				dimensions: settings.dimensions,
				points: (0..nb).map(|_| T::zero()).collect(),
			},
			fn_attraction: forces::choose_attraction(&settings),
			fn_gravity: forces::choose_gravity(&settings),
			fn_repulsion: Self::choose_repulsion(&settings),
			settings,
		}
	}

	#[inline]
	pub fn get_settings(&self) -> &Settings<T> {
		&self.settings
	}

	pub fn set_settings(&mut self, settings: Settings<T>) {
		self.fn_attraction = forces::choose_attraction(&settings);
		self.fn_gravity = forces::choose_gravity(&settings);
		self.fn_repulsion = Self::choose_repulsion(&settings);
		self.settings = settings;
	}

	/// Computes an iteration of ForceAtlas2
	pub fn iteration(&mut self) {
		self.init_iteration();
		self.apply_attraction();
		self.apply_repulsion();
		self.apply_gravity();
		self.apply_forces();
		self.apply_repulsion();
	}

	fn init_iteration(&mut self) {
		for (speed, old_speed) in self
			.speeds
			.points
			.iter_mut()
			.zip(self.old_speeds.points.iter_mut())
		{
			*old_speed = speed.clone();
			*speed = T::zero();
		}
	}

	#[inline]
	fn apply_attraction(&mut self) {
		(self.fn_attraction)(self)
	}

	#[inline]
	fn apply_gravity(&mut self) {
		(self.fn_gravity)(self)
	}

	#[inline]
	fn apply_repulsion(&mut self) {
		(self.fn_repulsion)(self)
	}

	fn apply_forces(&mut self) {
		for (pos, speed, old_speed) in izip!(
			self.points.iter_mut(),
			self.speeds.iter(),
			self.old_speeds.iter()
		) {
			let swinging = speed
				.iter()
				.zip(old_speed.iter())
				.map(|(s, old_s)| (s.clone() - old_s.clone()).pow_n(2u32))
				.sum::<T>();
			let traction = speed
				.iter()
				.zip(old_speed.iter())
				.map(|(s, old_s)| (s.clone() + old_s.clone()).pow_n(2u32))
				.sum::<T>();
			let f = (traction).ln() / (swinging.sqrt() + T::one());

			pos.iter_mut()
				.zip(speed.iter())
				.for_each(|(pos, speed)| *pos += speed.clone() * f.clone());
		}
	}
}

#[cfg(test)]
mod tests {
	use super::*;

	use alloc_counter::AllocCounterSystem;

	#[global_allocator]
	static A: AllocCounterSystem = AllocCounterSystem;

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
			vec![vec![-2.0, -2.0].as_slice(), vec![1.0, 2.0].as_slice()].into_iter(),
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
		assert_eq!(speed_1[0], 3.0);
		assert_eq!(speed_1[1], 4.0);
		assert_eq!(speed_2[0], -3.0);
		assert_eq!(speed_2[1], -4.0);

		layout.init_iteration();
		layout.apply_repulsion();

		let speed_1 = dbg!(layout.speeds.get(0));
		let speed_2 = dbg!(layout.speeds.get(1));

		assert!(speed_1[0] < 0.0);
		assert!(speed_1[1] < 0.0);
		assert!(speed_2[0] > 0.0);
		assert!(speed_2[1] > 0.0);
		assert_eq!(speed_1[0], -0.48);
		assert_eq!(speed_1[1], -0.64);
		assert_eq!(speed_2[0], 0.48);
		assert_eq!(speed_2[1], 0.64);

		layout.init_iteration();
		layout.apply_gravity();

		let speed_1 = dbg!(layout.speeds.get(0));
		let speed_2 = dbg!(layout.speeds.get(1));

		assert!(speed_1[0] > 0.0);
		assert!(speed_1[1] > 0.0);
		assert!(speed_2[0] < 0.0);
		assert!(speed_2[1] < 0.0);
		assert_eq!(speed_1[0], 2.0 / 2.0.sqrt());
		assert_eq!(speed_1[1], 2.0 / 2.0.sqrt());
		assert_eq!(speed_2[0], -2.0 / 5.0.sqrt());
		assert_eq!(speed_2[1], -4.0 / 5.0.sqrt());
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
				jitter_tolerance: 1.0,
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

	#[test]
	fn check_alloc() {
		let mut layout = Layout::<f64>::from_graph(
			vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 4), (3, 4)],
			5,
			Settings::default(),
		);

		layout.init_iteration();

		assert_eq!(
			alloc_counter::count_alloc(|| layout.apply_attraction()).0,
			(1, 0, 1)
		);
		assert_eq!(
			alloc_counter::count_alloc(|| layout.apply_repulsion()).0,
			(0, 0, 0)
		);
		assert_eq!(
			alloc_counter::count_alloc(|| layout.apply_gravity()).0,
			(0, 0, 0)
		);
	}
}
