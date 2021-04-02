use crate::util::*;

#[derive(Clone)]
pub struct Settings<T: Coord> {
	/// Number of spatial dimensions
	pub dimensions: usize,
	/// Move hubs (high degree nodes) to the center
	pub dissuade_hubs: bool,
	pub jitter_tolerance: T,
	/// Attraction coefficient
	pub ka: T,
	/// Gravity coefficient
	pub kg: T,
	/// Repulsion coefficient
	pub kr: T,
	/// Logarithmic attraction
	pub lin_log: bool,
	/// Prevent node overlapping for a prettier graph (node_size, kr_prime).
	///
	/// `node_size` is the radius around a node where the repulsion coefficient is `kr_prime`.
	/// `kr_prime` is arbitrarily set to `100.0` in Gephi implementation.
	pub prevent_overlapping: Option<(T, T)>,
	/// Gravity does not decrease with distance, resulting in a more compact graph.
	pub strong_gravity: bool,
	/// Optimize repulsion using Barnes-Hut algorithm (time passes from N^2 to NlogN).
	/// The argument is theta.
	///
	/// **Note**: only implemented for `T=f64` and `dimension` 2 or 3.
	#[cfg(feature = "barnes_hut")]
	pub barnes_hut: Option<T>,
}

impl<T: Coord> Default for Settings<T> {
	fn default() -> Self {
		Self {
			dimensions: 2,
			dissuade_hubs: false,
			jitter_tolerance: T::one(),
			ka: T::one(),
			kg: T::one(),
			kr: T::one(),
			lin_log: false,
			prevent_overlapping: None,
			strong_gravity: false,
			#[cfg(feature = "barnes_hut")]
			barnes_hut: None,
		}
	}
}

pub struct Layout<T: Coord> {
	pub edges: Vec<Edge>,
	pub nodes: Vec<Node>,
	/// List of the nodes' positions
	pub points: PointList<T>,
	pub(crate) settings: Settings<T>,
	pub speed: T,
	pub speed_efficiency: T,
	pub speeds: PointList<T>,
	pub old_speeds: PointList<T>,

	pub(crate) fn_attraction: fn(&mut Self),
	pub(crate) fn_gravity: fn(&mut Self),
	pub(crate) fn_repulsion: fn(&mut Self),
}
