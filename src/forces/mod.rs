pub mod attraction;
pub mod gravity;
pub mod repulsion;

use crate::{
	layout::{Layout, Settings},
	util::*,
};

#[doc(hidden)]
pub trait Repulsion<T: Coord + std::fmt::Debug> {
	fn choose_repulsion(settings: &Settings<T>) -> fn(&mut Layout<T>);
}

#[allow(clippy::collapsible_else_if)]
pub fn choose_attraction<T: Coord + std::fmt::Debug>(settings: &Settings<T>) -> fn(&mut Layout<T>) {
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
				attraction::apply_attraction
			}
		}
	}
}

pub fn choose_gravity<T: Coord + std::fmt::Debug>(settings: &Settings<T>) -> fn(&mut Layout<T>) {
	if settings.strong_gravity {
		gravity::apply_gravity_sg
	} else {
		gravity::apply_gravity
	}
}

default impl<T> Repulsion<T> for Layout<T>
//forces::Forces<T>
where
	T: Coord + std::fmt::Debug,
{
	fn choose_repulsion(settings: &Settings<T>) -> fn(&mut Layout<T>) {
		#[cfg(feature = "barnes_hut")]
		if settings.barnes_hut.is_some() {
			unimplemented!("Barnes-Hut only implemented for f64")
		}
		if settings.prevent_overlapping.is_some() {
			repulsion::apply_repulsion_po
		} else {
			repulsion::apply_repulsion
		}
	}
}

impl Repulsion<f64> for Layout<f64> /*forces::Forces<f64>*/ {
	fn choose_repulsion(settings: &Settings<f64>) -> fn(&mut Layout<f64>) {
		#[cfg(feature = "barnes_hut")]
		if settings.barnes_hut.is_some() {
			return match settings.dimensions {
				2 => {
					if settings.prevent_overlapping.is_some() {
						repulsion::apply_repulsion_bh_2d_po
					} else {
						repulsion::apply_repulsion_bh_2d
					}
				}
				3 => {
					if settings.prevent_overlapping.is_some() {
						repulsion::apply_repulsion_bh_3d_po
					} else {
						repulsion::apply_repulsion_bh_3d
					}
				}
				_ => unimplemented!("Barnes-Hut only implemented for 2D and 3D"),
			};
		}
		if settings.prevent_overlapping.is_some() {
			repulsion::apply_repulsion_po
		} else {
			match settings.dimensions {
				2 => repulsion::apply_repulsion_fast_2d,
				_ => repulsion::apply_repulsion_fast,
			}
		}
	}
}
