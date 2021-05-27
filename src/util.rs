use num_traits::cast::{ NumCast};
#[cfg(feature = "rand")]
use rand::Rng;
use num_traits::{Zero};
use std::ops::Neg;

pub type Coord = f32;

/// n-dimensional position
pub type Position = [Coord];

pub fn clone_slice_mut<T: Clone>(s: &[T]) -> Vec<T> {
	let mut v = valloc::<T>(s.len());
	let c: &mut [T] = v.as_mut_slice();
	/*for (i, e) in s.iter().enumerate() {
		c[i] = e.clone();
	}*/
	s.iter().zip(c.iter_mut()).for_each(|(e, c)| *c = e.clone());
	v
}

pub type Edge = (usize, usize);

pub enum Nodes {
	Mass(Vec<Coord>),
	Degree(usize),
}

pub fn norm(n: &Position) -> Coord {
	n.iter().map(|i| (*i).powi(2)).sum::<Coord>().sqrt()
}

/// Allocate Vec without initializing
pub fn valloc<T>(n: usize) -> Vec<T> {
	let mut v = Vec::with_capacity(n);
	unsafe {
		v.set_len(n);
	}
	v
}

pub(crate) unsafe fn split_at_mut_unchecked<T>(s: &mut [T], mid: usize) -> (&mut [T], &mut [T]) {
	let len = s.len();
	let ptr = s.as_mut_ptr();

	(
		std::slice::from_raw_parts_mut(ptr, mid),
		std::slice::from_raw_parts_mut(ptr.add(mid), len - mid),
	)
}

pub struct PointIter<'a> {
	pub dimensions: usize,
	pub offset: usize,
	pub list: &'a Vec<Coord>,
}

impl<'a> PointIter<'a> {
	/// Returns a raw pointer to the next element, and increments the counter by `n`.
	///
	/// # Safety
	/// Returned pointer may overflow the data.
	pub unsafe fn next_unchecked(&mut self, n: usize) -> *const Coord {
		let ptr = self.list.as_ptr().add(self.offset);
		self.offset += self.dimensions * n;
		ptr
	}
}

impl<'a> Iterator for PointIter<'a> {
	type Item = &'a [Coord];

	fn next(&mut self) -> Option<Self::Item> {
		if self.offset >= self.list.len() {
			return None;
		}
		let ret = unsafe {
			self.list
				.get_unchecked(self.offset..self.offset + self.dimensions)
		};
		self.offset += self.dimensions;
		Some(ret)
	}
}

pub struct PointIterMut<'a> {
	pub dimensions: usize,
	pub offset: usize,
	pub list: &'a mut Vec<Coord>,
}

impl<'a> Iterator for PointIterMut<'a> {
	type Item = &'a mut [Coord];

	fn next<'b>(&'b mut self) -> Option<Self::Item> {
		if self.offset >= self.list.len() {
			return None;
		}
		let ret: &'b mut [Coord] = unsafe {
			self.list
				.get_unchecked_mut(self.offset..self.offset + self.dimensions)
		};
		self.offset += self.dimensions;
		Some(unsafe { std::mem::transmute(ret) })
	}
}

#[derive(Clone)]
pub struct PointList {
	/// Number of coordinates in a vector
	pub dimensions: usize,
	/// List of the coordinates of the vectors
	pub points: Vec<Coord>,
}

impl<'a> PointList {
	pub fn get(&'a self, n: usize) -> &'a Position {
		let offset = n * self.dimensions;
		&self.points[offset..offset + self.dimensions]
	}

	/// # Safety
	/// `n` must be in bounds.
	pub unsafe fn get_unchecked(&'a self, n: usize) -> &'a Position {
		let offset = n * self.dimensions;
		&self.points.get_unchecked(offset..offset + self.dimensions)
	}

	pub fn get_clone(&self, n: usize) -> Vec<Coord> {
		clone_slice_mut(self.get(n))
	}

	pub fn get_clone_slice(&self, n: usize, v: &mut [Coord]) {
		v.clone_from_slice(self.get(n))
	}

	pub fn get_mut(&mut self, n: usize) -> &mut Position {
		let offset = n * self.dimensions;
		&mut self.points[offset..offset + self.dimensions]
	}

	/// n1 < n2
	pub fn get_2_mut(&mut self, n1: usize, n2: usize) -> (&mut Position, &mut Position) {
		let offset1 = n1 * self.dimensions;
		let offset2 = n2 * self.dimensions;
		unsafe {
			let (s1, s2) = split_at_mut_unchecked(&mut self.points, offset2);
			(
				s1.get_unchecked_mut(offset1..offset1 + self.dimensions),
				s2.get_unchecked_mut(..self.dimensions),
			)
		}
	}

	pub fn set(&mut self, n: usize, val: &Position) {
		let offset = n * self.dimensions;
		self.points[offset..offset + self.dimensions].clone_from_slice(val);
	}

	pub fn iter(&self) -> PointIter {
		PointIter {
			dimensions: self.dimensions,
			list: &self.points,
			offset: 0,
		}
	}

	pub fn iter_from(&self, offset: usize) -> PointIter {
		PointIter {
			dimensions: self.dimensions,
			list: &self.points,
			offset: offset * self.dimensions,
		}
	}

	pub fn iter_mut(&mut self) -> PointIterMut {
		PointIterMut {
			dimensions: self.dimensions,
			list: &mut self.points,
			offset: 0,
		}
	}

	pub fn iter_mut_from(&mut self, offset: usize) -> PointIterMut {
		PointIterMut {
			dimensions: self.dimensions,
			list: &mut self.points,
			offset: offset * self.dimensions,
		}
	}
}

/// Uniform random distribution of points on a n-sphere
///
/// `n` is the number of spatial dimensions (1 => two points; 2 => circle; 3 => sphere; etc.).
#[cfg(feature = "rand")]
pub fn _sample_unit_nsphere<R: Rng>(
	rng: &mut R,
	n: usize,
) -> Vec<Coord>
where
	rand::distributions::Standard: rand::distributions::Distribution<Coord>,
{
	let ray: Coord = NumCast::from(n).unwrap();
	let mut v = valloc(n);
	let mut d = Coord::zero();
	for x in v.iter_mut() {
		*x = rng.gen_range(ray.clone().neg()..ray.clone());
		d += (*x).powi(2);
	}
	d = d.sqrt();
	for x in v.iter_mut() {
		*x /= d.clone();
	}
	v
}

/// Uniform random distribution of points in a n-cube
///
/// `n` is the number of spatial dimensions (1 => segment; 2 => square; 3 => cube; etc.).
#[cfg(feature = "rand")]
pub fn sample_unit_ncube<R: Rng>(
	rng: &mut R,
	n: usize,
) -> Vec<Coord>
where
	rand::distributions::Standard: rand::distributions::Distribution<Coord>,
{
	let ray: Coord = NumCast::from(n).unwrap();
	let mut v = valloc(n);
	for x in v.iter_mut() {
		*x = rng.gen_range(ray.clone().neg()..ray.clone());
	}
	v
}

pub(crate) struct SendPtr<T>(pub std::ptr::NonNull<T>);

impl<T> Copy for SendPtr<T> {}
impl<T> Clone for SendPtr<T> {
	fn clone(&self) -> Self {
		SendPtr(self.0)
	}
}

unsafe impl<T> Send for SendPtr<T> {}
unsafe impl<T> Sync for SendPtr<T> {}

#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn test_clone_slice_mut() {
		let a = [1, 2, 3, 4, 5];
		let mut bv = clone_slice_mut(&a);
		let b = bv.as_mut_slice();
		b[2] = 6;
		assert_eq!(b.len(), 5);
		assert_eq!(b, [1, 2, 6, 4, 5]);
	}

	#[test]
	fn test_get_2_mut() {
		let mut a = PointList {
			dimensions: 2,
			points: vec![0., 1., 2., 3., 4., 5., 6., 7., 8., 9.],
		};
		let (s1, s2) = a.get_2_mut(1, 3);
		assert_eq!(s1.to_vec(), [2., 3.]);
		assert_eq!(s2.to_vec(), [6., 7.]);
	}
}
