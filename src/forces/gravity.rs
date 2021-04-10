use crate::{layout::Layout, util::*};

use itertools::izip;
use num_traits::cast::NumCast;

pub fn apply_gravity<T: Coord + std::fmt::Debug>(layout: &mut Layout<T>) {
	for (n, pos, speed) in izip!(
		layout.nodes.iter(),
		layout.points.iter(),
		layout.speeds.iter_mut()
	) {
		let d = norm(pos);
		if d.is_zero() {
			continue;
		}
		let f = <T as NumCast>::from(n.degree + 1).unwrap() * layout.settings.kg.clone() / d;
		for (speed, pos) in speed.iter_mut().zip(pos.iter()) {
			*speed -= f.clone() * pos.clone();
		}
	}
}

pub fn apply_gravity_sg<T: Coord + std::fmt::Debug>(layout: &mut Layout<T>) {
	for (n, pos, speed) in izip!(
		layout.nodes.iter(),
		layout.points.iter(),
		layout.speeds.iter_mut()
	) {
		let f = <T as NumCast>::from(n.degree + 1).unwrap() * layout.settings.kg.clone();
		for (speed, pos) in speed.iter_mut().zip(pos.iter()) {
			*speed -= f.clone() * pos.clone();
		}
	}
}
