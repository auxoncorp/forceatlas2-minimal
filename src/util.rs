use maths_traits::{
	algebra::group_like::{
		additive::Sub,
		multiplicative::{Div, DivAssign},
	},
	analysis::{ordered::Signed, RealExponential},
};
use num_traits::cast::{FromPrimitive, NumCast};
#[cfg(feature = "rand")]
use rand::Rng;

pub trait Coord = Clone
	+ Div<Self, Output = Self>
	+ DivAssign<Self>
	+ FromPrimitive
	+ NumCast
	+ From<f32>
	+ Signed
	+ RealExponential
	+ Sub<Self>
	+ std::iter::Sum;

/// n-dimensional position
pub type Position<T> = [T];

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

pub struct Node {
	pub degree: u32,
}

#[inline]
pub fn norm<T: Coord>(n: &Position<T>) -> T {
	n.iter().map(|i| i.clone().pow_n(2u32)).sum::<T>().sqrt()
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

impl<'a, T> PointIter<'a, T> {
	/// Returns a raw pointer to the next element, and increments the counter by `n`.
	///
	/// # Safety
	/// Returned pointer may overflow the data.
	pub unsafe fn next_unchecked(&mut self, n: usize) -> *const T {
		let ptr = self.list.as_ptr().add(self.offset);
		self.offset += self.dimensions * n;
		ptr
	}
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
	pub fn get_clone_slice(&self, n: usize, v: &mut [T]) {
		v.clone_from_slice(self.get(n))
	}
	#[inline(always)]
	pub fn get_mut(&mut self, n: usize) -> &mut Position<T> {
		let offset = n * self.dimensions;
		&mut self.points[offset..offset + self.dimensions]
	}
	/// n1 < n2
	#[inline(always)]
	pub fn get_2_mut(&mut self, n1: usize, n2: usize) -> (&mut Position<T>, &mut Position<T>) {
		let offset1 = n1 * self.dimensions;
		let offset2 = n2 * self.dimensions;
		let (s1, s2) = self.points.split_at_mut(offset2);
		(
			unsafe { s1.get_unchecked_mut(offset1..offset1 + self.dimensions) },
			unsafe { s2.get_unchecked_mut(..self.dimensions) },
		)
	}
	#[inline(always)]
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
pub fn _sample_unit_nsphere<T: Coord + Clone + DivAssign<T> + RealExponential, R: Rng>(
	rng: &mut R,
	n: usize,
) -> Vec<T>
where
	rand::distributions::Standard: rand::distributions::Distribution<T>,
	T: rand::distributions::uniform::SampleUniform + PartialOrd,
{
	let ray: T = NumCast::from(n).unwrap();
	let mut v = valloc(n);
	let mut d = T::zero();
	for x in v.iter_mut() {
		*x = rng.gen_range(ray.clone().neg()..ray.clone());
		d += x.clone().pow_n(2u32);
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
pub fn sample_unit_ncube<T: Coord + Clone + DivAssign<T> + RealExponential, R: Rng>(
	rng: &mut R,
	n: usize,
) -> Vec<T>
where
	rand::distributions::Standard: rand::distributions::Distribution<T>,
	T: rand::distributions::uniform::SampleUniform + PartialOrd,
{
	let ray: T = NumCast::from(n).unwrap();
	let mut v = valloc(n);
	for x in v.iter_mut() {
		*x = rng.gen_range(ray.clone().neg()..ray.clone());
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
