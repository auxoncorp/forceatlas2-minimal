pub mod attraction;
pub mod gravity;
pub mod repulsion;

use crate::{
	layout::{Layout, Settings},
};
use num_traits::Zero;

#[doc(hidden)]
pub trait Attraction {
	fn choose_attraction(settings: &Settings) -> fn(&mut Layout);
}

#[doc(hidden)]
pub trait Repulsion {
	fn choose_repulsion(settings: &Settings) -> fn(&mut Layout);
}

impl Attraction for Layout
{
	#[allow(clippy::collapsible_else_if)]
	fn choose_attraction(settings: &Settings) -> fn(&mut Layout) {
		if settings.prevent_overlapping.is_some() {
			if settings.lin_log {
				if settings.dissuade_hubs {
					attraction::apply_attraction_dh_log_po
				} else {
					attraction::apply_attraction_log_po
				}
			} else {
				if settings.dissuade_hubs {
					attraction::apply_attraction_dh_po
				} else {
					attraction::apply_attraction_po
				}
			}
		} else {
			if settings.lin_log {
				if settings.dissuade_hubs {
					attraction::apply_attraction_dh_log
				} else {
					attraction::apply_attraction_log
				}
			} else {
				if settings.dissuade_hubs {
					attraction::apply_attraction_dh
				} else {
					match settings.dimensions {
						2 => attraction::apply_attraction_2d,
						3 => attraction::apply_attraction_3d,
						_ => attraction::apply_attraction,
					}
				}
			}
		}
	}
}

pub fn choose_gravity(settings: &Settings) -> fn(&mut Layout) {
	if settings.kg.is_zero() {
		return |_| {};
	}
	if settings.strong_gravity {
		gravity::apply_gravity_sg
	} else {
		gravity::apply_gravity
	}
}


impl Repulsion for Layout {
	fn choose_repulsion(settings: &Settings) -> fn(&mut Layout) {
		if settings.prevent_overlapping.is_some() {
			repulsion::apply_repulsion_po
		} else {
			match settings.dimensions {
				2 => {
					#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
					{
						if is_x86_feature_detected!("avx2") {
							#[cfg(feature = "parallel")]
							if settings.chunk_size.is_some() {
								return repulsion::apply_repulsion_2d_simd_f32_parallel;
							}
							return repulsion::apply_repulsion_2d_simd_f32;
						}
					}
					#[cfg(feature = "parallel")]
					if settings.chunk_size.is_some() {
						return repulsion::apply_repulsion_2d_parallel;
					}
					repulsion::apply_repulsion_2d
				}
				3 => repulsion::apply_repulsion_3d,
				_ => {
					#[cfg(feature = "parallel")]
					if settings.chunk_size.is_some() {
						return repulsion::apply_repulsion_parallel;
					}
					repulsion::apply_repulsion
				}
			}
		}
	}
}
