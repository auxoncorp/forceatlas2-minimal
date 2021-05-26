use crate::{layout::Layout, util::*};

use std::marker::PhantomData;

#[cfg(feature = "parallel")]
pub use parallel::*;
#[cfg(all(feature = "parallel", any(target_arch = "x86", target_arch = "x86_64")))]
pub use parallel_simd::*;

pub struct Node<'a> {
	#[cfg(test)]
	pub ind: usize,
	pub mass: &'a Coord,
	pub n2_iter: NodeIter2<'a>,
	pub pos: &'a [Coord],
	pub speed: &'a mut [Coord],
}

pub struct NodeIter<'a> {
	pub ind: usize,
	pub(crate) layout: SendPtr<Layout>,
	pub offset: usize,
	pub(crate) _phantom: PhantomData<&'a mut Layout>,
}

impl <'a> NodeIter<'a> {
	unsafe fn unchecked_layout(&mut self) -> &'a mut Layout {
		core::mem::transmute(self.layout.0.as_mut())
	}
}

pub struct Node2<'a> {
	#[cfg(test)]
	pub ind: usize,
	pub mass: &'a Coord,
	pub pos: &'a [Coord],
	pub speed: &'a mut [Coord],
}

pub struct NodeIter2<'a> {
	pub ind: usize,
	pub(crate) layout: SendPtr<Layout>,
	pub offset: usize,
	pub(crate) _phantom: PhantomData<&'a mut Layout>,
}

impl <'a> NodeIter2<'a> {
	unsafe fn unchecked_layout(&mut self) -> &'a mut Layout {
		core::mem::transmute(self.layout.0.as_mut())
	}
}

impl<'a> Iterator for NodeIter<'a> {
	type Item = Node<'a>;

	fn next(&mut self) -> Option<Self::Item> {
		let layout: &'a Layout = unsafe { self.unchecked_layout() };
		if let Some(mass) = layout.masses.get(self.ind) {
			Some({
				let next_offset = self.offset + layout.settings.dimensions;
				let next_ind = self.ind + 1;
				let ret = Node {
					#[cfg(test)]
					ind: self.ind,
					mass,
					n2_iter: NodeIter2 {
						ind: next_ind,
						layout: self.layout,
						offset: next_offset,
						_phantom: PhantomData,
					},
					pos: unsafe { layout.points.points.get_unchecked(self.offset..next_offset) },
					speed: unsafe {
						self.unchecked_layout()
							.speeds
							.points
							.get_unchecked_mut(self.offset..next_offset)
					},
				};
				self.offset = next_offset;
				self.ind = next_ind;
				ret
			})
		} else {
			None
		}
	}
}

impl<'a> Iterator for NodeIter2<'a> {
	type Item = Node2<'a>;

	fn next(&mut self) -> Option<Self::Item> {
		let layout = unsafe { self.unchecked_layout() };
		if let Some(mass) = layout.masses.get(self.ind) {
			Some({
				let next_offset = self.offset + layout.settings.dimensions;
				let ret = Node2 {
					#[cfg(test)]
					ind: self.offset / layout.settings.dimensions,
					mass,
					pos: unsafe { layout.points.points.get_unchecked(self.offset..next_offset) },
					speed: unsafe {
						self.unchecked_layout()
							.speeds
							.points
							.get_unchecked_mut(self.offset..next_offset)
					},
				};
				self.offset = next_offset;
				self.ind += 1;
				ret
			})
		} else {
			None
		}
	}
}

#[cfg(feature = "parallel")]
mod parallel {
	use super::*;

	pub struct NodePar<'a> {
		#[cfg(test)]
		pub ind: usize,
		pub mass: &'a Coord,
		pub n2_iter: NodeParIter2<'a>,
		pub pos: &'a [Coord],
		pub speed: &'a mut [Coord],
	}

	pub struct NodeParIter<'a> {
		pub end: usize,
		pub ind: usize,
		pub(crate) layout: SendPtr<Layout>,
		pub n2_start: usize,
		pub n2_start_ind: usize,
		pub n2_end: usize,
		pub offset: usize,
		pub(crate) _phantom: PhantomData<&'a mut Layout>,
	}
	impl <'a> NodeParIter<'a> {
		unsafe fn unchecked_layout(&mut self) -> &'a mut Layout {
			core::mem::transmute(self.layout.0.as_mut())
		}
	}

	pub struct NodePar2<'a> {
		#[cfg(test)]
		pub ind: usize,
		pub mass: &'a Coord,
		pub pos: &'a [Coord],
		pub speed: &'a mut [Coord],
	}

	pub struct NodeParIter2<'a> {
		pub end: usize,
		pub ind: usize,
		pub(crate) layout: SendPtr<Layout>,
		pub offset: usize,
		pub(crate) _phantom: PhantomData<&'a mut Layout>,
	}
	impl <'a> NodeParIter2<'a> {
		unsafe fn unchecked_layout(&mut self) -> &'a mut Layout {
			core::mem::transmute(self.layout.0.as_mut())
		}
	}

	impl<'a> Iterator for NodeParIter<'a> {
		type Item = NodePar<'a>;

		fn next(&mut self) -> Option<Self::Item> {
			if self.offset < self.end {
				Some({
					let layout = unsafe { self.unchecked_layout() };
					let next_offset = self.offset + layout.settings.dimensions;
					let next_ind = self.ind + 1;
					let ret = NodePar {
						#[cfg(test)]
						ind: self.ind,
						mass: unsafe { layout.masses.get_unchecked(self.ind) },
						n2_iter: NodeParIter2 {
							end: self.n2_end,
							ind: self.n2_start_ind.max(next_ind),
							layout: self.layout,
							offset: self.n2_start.max(next_offset),
							_phantom: PhantomData::default(),
						},
						pos: unsafe {
							layout.points.points.get_unchecked(self.offset..next_offset)
						},
						speed: unsafe {
							self.unchecked_layout()
								.speeds
								.points
								.get_unchecked_mut(self.offset..next_offset)
						},
					};
					self.offset = next_offset;
					self.ind = next_ind;
					ret
				})
			} else {
				None
			}
		}
	}

	impl<'a> Iterator for NodeParIter2<'a> {
		type Item = NodePar2<'a>;

		fn next(&mut self) -> Option<Self::Item> {
			if self.offset < self.end {
				Some({
					let layout = unsafe { self.unchecked_layout() };
					let next_offset = self.offset + layout.settings.dimensions;
					let ret = NodePar2 {
						#[cfg(test)]
						ind: self.ind,
						mass: unsafe { layout.masses.get_unchecked(self.ind) },
						pos: unsafe {
							layout.points.points.get_unchecked(self.offset..next_offset)
						},
						speed: unsafe {
							self.unchecked_layout()
								.speeds
								.points
								.get_unchecked_mut(self.offset..next_offset)
						},
					};
					self.offset = next_offset;
					self.ind += 1;
					ret
				})
			} else {
				None
			}
		}
	}
}

#[cfg(all(feature = "parallel", any(target_arch = "x86", target_arch = "x86_64")))]
mod parallel_simd {
	use super::*;

	pub struct NodeParSimd<'a, const N: usize> {
		pub ind: usize,
		pub mass: &'a Coord,
		pub n2_iter: NodeParSimdIter2<'a, N>,
		pub pos: &'a [Coord],
		pub speed: *mut Coord,
	}

	pub struct NodeParSimdIter<'a, const N: usize> {
		pub end: usize,
		pub ind: usize,
		pub(crate) layout: SendPtr<Layout>,
		pub n2_start: usize,
		pub n2_start_ind: usize,
		pub n2_end: usize,
		pub offset: usize,
		pub(crate) _phantom: PhantomData<&'a mut Layout>,
	}
	impl <'a, const N: usize> NodeParSimdIter<'a, N> {
		unsafe fn unchecked_layout(&mut self) -> &'a Layout {
			core::mem::transmute(self.layout.0.as_mut())
		}
	}

	pub struct NodeParSimd2<const N: usize> {
		#[cfg(test)]
		pub ind: usize,
		pub mass: *mut Coord,
		pub pos: *mut Coord,
		pub speed: *mut Coord,
	}

	pub struct NodeParSimdIter2<'a, const N: usize> {
		pub end: usize,
		pub ind: usize,
		pub(crate) layout: SendPtr<Layout>,
		pub offset: usize,
		pub(crate) _phantom: PhantomData<&'a mut Layout>,
	}

	impl<'a, const N: usize> Iterator for NodeParSimdIter<'a, N> {
		type Item = NodeParSimd<'a, N>;

		fn next(&mut self) -> Option<Self::Item> {
			if self.offset < self.end {
				Some({
					let layout = unsafe { self.unchecked_layout() };
					let next_offset = self.offset + layout.settings.dimensions;
					let next_ind = self.ind + 1;
					let n2_start_ind = self.n2_start_ind.max(next_ind);
					let ret = NodeParSimd {
						ind: self.ind,
						mass: unsafe { layout.masses.get_unchecked(self.ind) },
						n2_iter: NodeParSimdIter2 {
							end: self.n2_end - (self.n2_end - n2_start_ind) % N,
							ind: n2_start_ind,
							layout: self.layout,
							offset: self.n2_start.max(next_offset),
							_phantom: PhantomData::default(),
						},
						pos: unsafe {
							layout.points.points.get_unchecked(self.offset..next_offset)
						},
						speed: unsafe {
							self.layout
								.0
								.as_mut()
								.speeds
								.points
								.as_mut_ptr()
								.add(self.offset)
						},
					};
					self.offset = next_offset;
					self.ind = next_ind;
					ret
				})
			} else {
				None
			}
		}
	}

	impl<'a, const N: usize> Iterator for &mut NodeParSimdIter2<'a, N> {
		type Item = NodeParSimd2<N>;

		fn next(&mut self) -> Option<Self::Item> {
			if self.offset < self.end {
				Some({
					let layout = unsafe { self.layout.0.as_mut() };
					let next_offset = self.offset + N * layout.settings.dimensions;
					let ret = NodeParSimd2 {
						#[cfg(test)]
						ind: self.ind,
						mass: unsafe { layout.masses.as_mut_ptr().add(self.ind) },
						pos: unsafe { layout.points.points.as_mut_ptr().add(self.offset) },
						speed: unsafe {
							self.layout
								.0
								.as_mut()
								.speeds
								.points
								.as_mut_ptr()
								.add(self.offset)
						},
					};
					self.offset = next_offset;
					self.ind += N;
					ret
				})
			} else {
				None
			}
		}
	}
}
