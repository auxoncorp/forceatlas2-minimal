use maths_traits::{
	algebra::group_like::{
		additive::Sub,
		multiplicative::{Div, DivAssign},
	},
	analysis::{ordered::Signed, RealExponential},
};
#[cfg(feature = "rand")]
use rand::Rng;

pub trait Coord = Clone
	+ Div<Self, Output = Self>
	+ DivAssign<Self>
	+ From<u32>
	+ Signed
	+ RealExponential
	+ Sub<Self>;

/// n-dimensional position
pub type Position<T> = [T];

pub fn clone_slice_mut<'a, T: Clone>(s: &'a [T]) -> Vec<T> {
	let mut v = valloc::<T>(s.len());
	let c: &mut [T] = v.as_mut_slice();
	for (i, e) in s.iter().enumerate() {
		c[i] = e.clone();
	}
	v
}

pub type Edge = (usize, usize);

pub struct Node {
	pub degree: u32,
}

#[inline]
pub fn norm<T: Coord>(n: &Position<T>) -> T {
	let mut sum = T::zero();
	for i in n.iter() {
		sum += i.clone().pow_n(2u32);
	}
	sum.sqrt()
}

/// Allocate Vec without initializing
#[inline]
pub fn valloc<T>(n: usize) -> Vec<T> {
	let mut v = Vec::with_capacity(n);
	unsafe {
		v.set_len(n);
	}
	v
}

pub struct PointIter<'a, T> {
	pub dimensions: usize,
	pub offset: usize,
	pub list: &'a Vec<T>,
}

impl<'a, T> Iterator for PointIter<'a, T> {
	type Item = &'a [T];

	fn next(&mut self) -> Option<Self::Item> {
		if self.offset >= self.list.len() {
			return None;
		}
		let ret = &self.list[self.offset..self.offset + self.dimensions];
		self.offset += self.dimensions;
		Some(ret)
	}
}

pub struct PointIterMut<'a, T> {
	pub dimensions: usize,
	pub offset: usize,
	pub list: &'a mut Vec<T>,
}

impl<'a, T> Iterator for PointIterMut<'a, T> {
	type Item = &'a mut [T];

	fn next<'b>(&'b mut self) -> Option<Self::Item> {
		if self.offset >= self.list.len() {
			return None;
		}
		let ret: &'b mut [T] = &mut self.list[self.offset..self.offset + self.dimensions];
		self.offset += self.dimensions;
		Some(unsafe { std::mem::transmute(ret) })
	}
}

pub struct PointList<T: Coord> {
	pub dimensions: usize,
	pub points: Vec<T>,
}

impl<'a, T: Coord> PointList<T> {
	#[inline]
	pub fn get(&'a self, n: usize) -> &'a Position<T> {
		let offset = n * self.dimensions;
		&self.points[offset..offset + self.dimensions]
	}
	#[inline]
	pub fn get_clone(&self, n: usize) -> Vec<T> {
		clone_slice_mut(self.get(n))
	}
	#[inline]
	pub fn get_mut(&mut self, n: usize) -> &mut Position<T> {
		let offset = n * self.dimensions;
		&mut self.points[offset..offset + self.dimensions]
	}
	#[inline]
	pub fn set(&mut self, n: usize, val: &Position<T>) {
		let offset = n * self.dimensions;
		self.points[offset..offset + self.dimensions].clone_from_slice(val);
	}
	pub fn iter(&self) -> PointIter<T> {
		PointIter {
			dimensions: self.dimensions,
			list: &self.points,
			offset: 0,
		}
	}
	pub fn iter_mut(&mut self) -> PointIterMut<T> {
		PointIterMut {
			dimensions: self.dimensions,
			list: &mut self.points,
			offset: 0,
		}
	}
}

/// Uniform random distribution of points on a n-sphere
///
/// `n` is the number of spatial dimensions (1 => two points; 2 => circle; 3 => sphere; etc.).
#[cfg(feature = "rand")]
pub fn sample_unit_nsphere<T: Clone + DivAssign<T> + RealExponential, R: Rng>(
	rng: &mut R,
	n: usize,
) -> Vec<T>
where
	rand::distributions::Standard: rand::distributions::Distribution<T>,
	T: rand::distributions::uniform::SampleUniform,
{
	let mut v = valloc(n);
	let mut d = T::zero();
	for x in v.iter_mut() {
		*x = rng.gen_range(T::one().neg(), T::one());
		d += x.clone().pow_n(2u32);
	}
	d = d.sqrt();
	for x in v.iter_mut() {
		*x /= d.clone();
	}
	v
}

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
}
