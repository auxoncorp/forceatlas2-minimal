use crate::{iter::*, util::*};

use rayon::prelude::*;
use std::marker::PhantomData;
use num_traits::One;

#[derive(Clone)]
pub struct Settings {
	/// Optimize repulsion using Barnes-Hut algorithm (time passes from N^2 to NlogN)
	/// The argument is theta.
	///
	/// **Note**: only implemented for `T=f64` and `dimension` 2 or 3.
	#[cfg(feature = "barnes_hut")]
	pub barnes_hut: Option<Coord>,
	/// Number of nodes computed by each thread
	///
	/// Only used in repulsion computation. Set to `None` to turn off parallelization.
	/// This number should be big enough to minimize thread management,
	/// but small enough to maximize concurrency.
	///
	/// Requires `T: Send + Sync`
	#[cfg(feature = "parallel")]
	pub chunk_size: Option<usize>,
	/// Number of spatial dimensions
	pub dimensions: usize,
	/// Move hubs (high degree nodes) to the center
	pub dissuade_hubs: bool,
	/// Attraction coefficient
	pub ka: Coord,
	/// Gravity coefficient
	pub kg: Coord,
	/// Repulsion coefficient
	pub kr: Coord,
	/// Logarithmic attraction
	pub lin_log: bool,
	/// Prevent node overlapping for a prettier graph (node_size, kr_prime).
	///
	/// `node_size` is the radius around a node where the repulsion coefficient is `kr_prime`.
	/// `kr_prime` is arbitrarily set to `100.0` in Gephi implementation.
	pub prevent_overlapping: Option<(Coord, Coord)>,
	/// Gravity does not decrease with distance, resulting in a more compact graph.
	pub strong_gravity: bool,
}

impl Default for Settings {
	fn default() -> Self {
		Self {
			#[cfg(feature = "barnes_hut")]
			barnes_hut: None,
			#[cfg(feature = "parallel")]
			chunk_size: Some(256),
			dimensions: 2,
			dissuade_hubs: false,
			ka: Coord::one(),
			kg: Coord::one(),
			kr: Coord::one(),
			lin_log: false,
			prevent_overlapping: None,
			strong_gravity: false,
		}
	}
}

pub struct Layout {
	pub edges: Vec<Edge>,
	pub masses: Vec<Coord>,
	/// List of the nodes' positions
	pub points: PointList,
	pub(crate) settings: Settings,
	pub speeds: PointList,
	pub old_speeds: PointList,

	pub(crate) fn_attraction: fn(&mut Self),
	pub(crate) fn_gravity: fn(&mut Self),
	pub(crate) fn_repulsion: fn(&mut Self),
}

impl Layout {
	pub fn iter_nodes(&mut self) -> NodeIter {
		NodeIter {
			ind: 0,
			layout: SendPtr(self.into()),
			offset: 0,
			_phantom: PhantomData::default(),
		}
	}
}

#[cfg(feature = "parallel")]
impl Layout {
	pub fn iter_par_nodes(
		&mut self,
		chunk_size: usize,
	) -> impl Iterator<Item = impl ParallelIterator<Item = NodeParIter>> {
		let ptr = SendPtr(self.into());
		let dimensions = self.settings.dimensions;
		let chunk_size_d = chunk_size * dimensions;
		let n = self.masses.len() * dimensions;
		(0..self.masses.len()).step_by(chunk_size).map(move |y0| {
			let y0_d = y0 * dimensions;
			(0..self.masses.len() - y0)
				.into_par_iter()
				.step_by(chunk_size)
				.map(move |x0| {
					let x0_d = x0 * dimensions;
					NodeParIter {
						end: (x0_d + chunk_size_d).min(n),
						ind: x0,
						layout: ptr,
						n2_start: x0_d + y0_d,
						n2_start_ind: x0 + y0,
						n2_end: (x0_d + y0_d + chunk_size_d).min(n),
						offset: x0_d,
						_phantom: PhantomData::default(),
					}
				})
		})
	}
}

#[cfg(all(feature = "parallel", any(target_arch = "x86", target_arch = "x86_64")))]
impl Layout {
	pub fn iter_par_simd_nodes<const N: usize>(
		&mut self,
		chunk_size: usize,
	) -> impl Iterator<Item = impl ParallelIterator<Item = NodeParSimdIter<N>>> {
		let ptr = SendPtr(self.into());
		let dimensions = self.settings.dimensions;
		let chunk_size_d = chunk_size * dimensions;
		let n = self.masses.len() * dimensions;
		(0..self.masses.len()).step_by(chunk_size).map(move |y0| {
			let y0_d = y0 * dimensions;
			(0..self.masses.len() - y0)
				.into_par_iter()
				.step_by(chunk_size)
				.map(move |x0| {
					let x0_d = x0 * dimensions;
					NodeParSimdIter {
						end: (x0_d + chunk_size_d).min(n),
						ind: x0,
						layout: ptr,
						n2_start: x0_d + y0_d,
						n2_start_ind: x0 + y0,
						n2_end: (x0_d + y0_d + chunk_size_d).min(n),
						offset: x0_d,
						_phantom: PhantomData::default(),
					}
				})
		})
	}
}

#[cfg(test)]
mod test {
	use super::*;
	use itertools::iproduct;
	use std::collections::BTreeSet;
	#[cfg(feature = "parallel")]
	use std::sync::{Arc, RwLock};

	#[test]
	fn test_iter_nodes() {
		for n_nodes in 1usize..16 {
			let mut layout =
				Layout::from_graph(vec![], Nodes::Degree(n_nodes), Settings::default());
			let mut hits = iproduct!(0..n_nodes, 0..n_nodes)
				.filter(|(n1, n2)| n1 < n2)
				.collect::<BTreeSet<(usize, usize)>>();
			let points = layout.points.clone();
			for n1 in layout.iter_nodes() {
				for n2 in n1.n2_iter {
					assert!(hits.remove(&(n1.ind, n2.ind)));
					assert_eq!(n1.pos, points.get(n1.ind));
					assert_eq!(n2.pos, points.get(n2.ind));
				}
			}
			assert!(hits.is_empty());
		}
	}

	#[test]
	#[cfg(feature = "parallel")]
	fn test_iter_par_nodes() {
		for n_nodes in 1usize..16 {
			let mut layout = Layout::from_graph(
				vec![],
				Nodes::Mass((1..n_nodes + 1).map(|i| i as f32).collect()),
				Settings::default(),
			);
			layout
				.speeds
				.iter_mut()
				.enumerate()
				.for_each(|(i, speed)| speed.iter_mut().for_each(|speed| *speed = i as f32));
			let hits = Arc::new(RwLock::new(
				iproduct!(0..n_nodes, 0..n_nodes)
					.filter(|(n1, n2)| n1 < n2)
					.collect::<BTreeSet<(usize, usize)>>(),
			));
			let points = layout.points.clone();
			let speeds = layout.speeds.clone();
			for chunk_iter in layout.iter_par_nodes(4) {
				chunk_iter.for_each(|n1_iter| {
					for n1 in n1_iter {
						for n2 in n1.n2_iter {
							let mut hits = hits.write().unwrap();
							assert!(hits.remove(&(n1.ind, n2.ind)));
							assert_eq!(n1.pos, points.get(n1.ind));
							assert_eq!(n2.pos, points.get(n2.ind));
							assert_eq!(n1.speed, speeds.get(n1.ind));
							assert_eq!(n2.speed, speeds.get(n2.ind));
							assert_eq!(*n1.mass, n1.ind as f32 + 1.);
							assert_eq!(*n2.mass, n2.ind as f32 + 1.);
						}
					}
				});
			}
			assert!(hits.read().unwrap().is_empty());
		}
	}

	#[test]
	#[cfg(feature = "parallel")]
	fn test_iter_par_simd_nodes() {
		for n_nodes in 1usize..32 {
			let mut layout = Layout::from_graph(
				vec![],
				Nodes::Mass((1..n_nodes + 1).map(|i| i as f32).collect()),
				Settings::default(),
			);
			layout
				.speeds
				.iter_mut()
				.enumerate()
				.for_each(|(i, speed)| speed.iter_mut().for_each(|speed| *speed = i as f32));
			let hits = Arc::new(RwLock::new(
				iproduct!(0..n_nodes, 0..n_nodes)
					.filter(|(n1, n2)| n1 < n2)
					.collect::<BTreeSet<(usize, usize)>>(),
			));
			let points = layout.points.clone();
			let speeds = layout.speeds.clone();
			for chunk_iter in layout.iter_par_simd_nodes::<4>(16) {
				chunk_iter.for_each(|n1_iter| {
					for mut n1 in n1_iter {
						for n2 in &mut n1.n2_iter {
							let mut hits = hits.write().unwrap();
							assert!(hits.remove(&(n1.ind, n2.ind)));
							assert!(hits.remove(&(n1.ind, n2.ind + 1)));
							assert!(hits.remove(&(n1.ind, n2.ind + 2)));
							assert!(hits.remove(&(n1.ind, n2.ind + 3)));
							assert_eq!(n1.pos, points.get(n1.ind));
							assert_eq!(n2.pos as *const f32, points.get(n2.ind).as_ptr());
							assert_eq!(n2.pos as *const f32, unsafe {
								points.get(n2.ind).as_ptr().add(1)
							});
							assert_eq!(n2.pos as *const f32, unsafe {
								points.get(n2.ind).as_ptr().add(2)
							});
							assert_eq!(n2.pos as *const f32, unsafe {
								points.get(n2.ind).as_ptr().add(3)
							});
							assert_eq!(n1.speed as *const f32, speeds.get(n1.ind).as_ptr());
							assert_eq!(n2.speed as *const f32, speeds.get(n2.ind).as_ptr());
							assert_eq!(n2.speed as *const f32, unsafe {
								speeds.get(n2.ind).as_ptr().add(1)
							});
							assert_eq!(n2.speed as *const f32, unsafe {
								speeds.get(n2.ind).as_ptr().add(2)
							});
							assert_eq!(n2.speed as *const f32, unsafe {
								speeds.get(n2.ind).as_ptr().add(3)
							});
							assert_eq!(*n1.mass, n1.ind as f32 + 1.);
							assert_eq!(unsafe { *n2.mass }, n2.ind as f32 + 1.);
							assert_eq!(unsafe { *n2.mass.add(1) }, n2.ind as f32 + 1.);
							assert_eq!(unsafe { *n2.mass }, n2.ind as f32 + 1.);
							assert_eq!(unsafe { *n2.mass }, n2.ind as f32 + 1.);
						}

						for n2 in n1.n2_iter.ind
							..(n1.n2_iter.ind + (n1.n2_iter.ind - n1.ind + 1) / 4 * 4 + 4)
								.min(n_nodes)
						{
							let mut hits = hits.write().unwrap();
							println!("rem {} {}", n1.ind, n2);
							assert!(hits.remove(&(n1.ind, n2)));
						}
					}
				});
			}
			println!("{}", n_nodes);
			println!("{:?}", hits);
			assert!(hits.read().unwrap().is_empty());
		}
	}
}
