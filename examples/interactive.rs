use forceatlas2::*;
use itertools::izip;
use plotters::prelude::*;
use sdl2;
use sdl2::{
	event::{Event, WindowEvent},
	keyboard::Keycode,
	pixels::PixelFormatEnum,
};
use std::{
	io::BufRead,
	sync::{Arc, RwLock},
	thread,
	time::Duration,
};
use std::ops::AddAssign;

const DRAW_LINKS: bool = true;
const SIZE: (u32, u32) = (800, 800);
const FRAMEDUR: u64 = 30;
type T = f32;

pub fn run(image: Arc<RwLock<(u32, u32, Vec<u8>)>>, size: Arc<RwLock<(u32, u32)>>) {
	let sdl_context = sdl2::init().unwrap();
	let video_subsystem = sdl_context.video().unwrap();

	let window = {
		let size = size.read().unwrap();
		video_subsystem
			.window("ForceAtlas2", size.0, size.1)
			.resizable()
			.position_centered()
			.opengl()
			.build()
			.unwrap()
	};

	let mut canvas = window.into_canvas().build().unwrap();
	let texture_creator = canvas.texture_creator();

	let mut event_pump = sdl_context.event_pump().unwrap();
	let interval = Duration::from_millis(FRAMEDUR);

	'running: loop {
		for event in event_pump.poll_iter() {
			match event {
				Event::Quit { .. }
				| Event::KeyDown {
					keycode: Some(Keycode::Escape),
					..
				} => {
					break 'running;
				}
				Event::Window {
					win_event: WindowEvent::SizeChanged(w, h),
					..
				} => {
					let mut size = size.write().unwrap();
					*size = (w as u32, h as u32);
				}
				_ => {}
			}
		}

		let image = image.read().unwrap();
		let mut texture = texture_creator
			.create_texture_streaming(PixelFormatEnum::RGB24, image.0, image.1)
			.unwrap();
		texture
			.update(None, &image.2, (image.0 * 3u32) as usize)
			.unwrap();
		canvas.copy(&texture, None, None).unwrap();

		canvas.present();

		std::thread::sleep(interval);
	}
}

fn main() {
	let file = std::fs::File::open(
		std::env::args()
			.skip(1)
			.next()
			.expect("Usage: csv_import <csv_file>"),
	)
	.expect("Cannot open file");

	let iters = Arc::new(RwLock::new(0u32));
	let mut nodes = 0usize;
	let mut edges = Vec::<(usize, usize)>::new();
	for (i, line) in std::io::BufReader::new(file).lines().enumerate() {
		let line = line.expect("Error reading CSV");
		let mut columns = line.split(&[' ', '\t', ',', ';'][..]);
		if let (Some(n1), Some(n2)) = (columns.next(), columns.skip_while(|&c| c.is_empty()).next())
		{
			if let (Ok(n1), Ok(n2)) = (n1.parse(), n2.parse()) {
				if n1 > nodes {
					nodes = n1;
				}
				if n2 > nodes {
					nodes = n2;
				}
				if n1 != n2 {
					edges.push(if n1 < n2 { (n1, n2) } else { (n2, n1) });
				}
			} else {
				eprintln!("Ignored line {} has bad number format", i);
			}
		} else {
			eprintln!("Ignored line {} has <2 columns", i);
		}
	}
	nodes += 1;

	println!("Nodes: {}", nodes);

	let mut settings = Settings {
		chunk_size: Some(256),
		dimensions: 2,
		dissuade_hubs: false,
		ka: 0.01,
		kg: 0.001,
		kr: 0.002,
		lin_log: false,
		prevent_overlapping: None,
		strong_gravity: false,
	};

	let layout = Arc::new(RwLock::new(Layout::from_graph(
		edges.clone(),
		Nodes::Degree(nodes),
		settings.clone(),
	)));
	let size = Arc::new(RwLock::new(SIZE));
	let image = Arc::new(RwLock::new((
		SIZE.0,
		SIZE.1,
		vec![0u8; SIZE.0 as usize * SIZE.1 as usize * 3],
	)));
	let computing = Arc::new(RwLock::new(false));
	let sleep = Arc::new(RwLock::new(50u64));
	draw_graph(layout.clone(), image.clone(), size.clone());

	thread::spawn({
		let image = image.clone();
		let computing = computing.clone();
		let layout = layout.clone();
		let sleep = sleep.clone();
		let size = size.clone();

		let interval = Duration::from_millis(FRAMEDUR);
		move || loop {
			if *computing.read().unwrap() {
				draw_graph(layout.clone(), image.clone(), size.clone());
				thread::sleep(interval);
			} else {
				thread::sleep(std::time::Duration::from_millis(*sleep.read().unwrap()));
			}
		}
	});

	thread::spawn({
		let computing = computing.clone();
		let layout = layout.clone();
		let iters = iters.clone();
		move || loop {
			if *computing.read().unwrap() {
				layout.write().unwrap().iteration();
				iters.write().unwrap().add_assign(1);
			} else {
				thread::sleep(std::time::Duration::from_millis(*sleep.read().unwrap()));
			}
		}
	});

	thread::spawn(move || run(image, size));

	thread::spawn(move || loop {
		let cmd = scanrs::scanln();
		let mut args = cmd.split(' ');
		match args.next() {
			Some("q") => break,
			Some("h") => {
				println!("h          help");
				println!("i          display number of iterations");
				println!("p          display parameters");
				println!("q          quit");
				println!("r          reset");
				println!("s          start/stop");
				println!("cs [usize] chunk size (leave empty to disable parallel)");
				println!("ka <f64>   attraction");
				println!("kg <f64>   gravity");
				println!("kr <f64>   repulsion");
			}
			Some("s") => {
				let mut computing = computing.write().unwrap();
				*computing = !*computing;
			}
			Some("cs") => {
				settings.chunk_size = args.next().map(|chunk_size| chunk_size.parse().unwrap());
				layout.write().unwrap().set_settings(settings.clone());
			}
			Some("ka") => {
				settings.ka = args.next().unwrap().parse().unwrap();
				layout.write().unwrap().set_settings(settings.clone());
			}
			Some("kg") => {
				settings.kg = args.next().unwrap().parse().unwrap();
				layout.write().unwrap().set_settings(settings.clone());
			}
			Some("kr") => {
				settings.kr = args.next().unwrap().parse().unwrap();
				layout.write().unwrap().set_settings(settings.clone());
			}
			Some("r") => {
				let mut layout = layout.write().unwrap();
				let mut iters = iters.write().unwrap();
				*iters = 0;
				*layout = Layout::from_graph(edges.clone(), Nodes::Degree(nodes), settings.clone());
			}
			Some("p") => println!(
				"ka={}  kg={}  kr={}  cs={:?}",
				settings.ka, settings.kg, settings.kr, settings.chunk_size
			),
			Some("i") => println!("{}", iters.read().unwrap()),
			_ => println!("Unknown command"),
		}
	})
	.join()
	.unwrap();
}

fn draw_graph(
	layout: Arc<RwLock<Layout>>,
	image: Arc<RwLock<(u32, u32, Vec<u8>)>>,
	size: Arc<RwLock<(u32, u32)>>,
) {
	let mut image = image.write().unwrap();
	{
		let size = size.read().unwrap();
		if size.0 * size.1 != image.0 * image.1 {
			image.2.resize((size.0 * size.1 * 3) as usize, 0);
		}
		image.0 = size.0;
		image.1 = size.1;
	}
	let size = (image.0, image.1);
	let root = BitMapBackend::with_buffer(&mut image.2, size).into_drawing_area();
	root.fill(&WHITE).unwrap();

	let layout = layout.read().unwrap();

	let mut min_v = layout.points.get_clone(0);
	let mut max_v = min_v.clone();
	let min = min_v.as_mut_slice();
	let max = max_v.as_mut_slice();
	for pos in layout.points.iter() {
		if pos[0] < min[0] {
			min[0] = pos[0];
		}
		if pos[1] < min[1] {
			min[1] = pos[1];
		}
		if pos[0] > max[0] {
			max[0] = pos[0];
		}
		if pos[1] > max[1] {
			max[1] = pos[1];
		}
	}
	let graph_size = (max[0] - min[0], max[1] - min[1]);
	let factor = {
		let factors = (size.0 as T / graph_size.0, size.1 as T / graph_size.1);
		if factors.0 > factors.1 {
			min[0] -= (size.0 as T / factors.1 - graph_size.0) / 2.0;
			factors.1
		} else {
			min[1] -= (size.1 as T / factors.0 - graph_size.1) / 2.0;
			factors.0
		}
	};

	if DRAW_LINKS {
		let link_color = RGBColor(5, 5, 5).mix(0.05);

		for (h1, h2) in layout.edges.iter() {
			root.draw(&PathElement::new(
				vec![
					{
						let pos = layout.points.get(*h1);
						unsafe {
							(
								((pos[0] - min[0]) * factor).to_int_unchecked::<i32>(),
								((pos[1] - min[1]) * factor).to_int_unchecked::<i32>(),
							)
						}
					},
					{
						let pos = layout.points.get(*h2);
						unsafe {
							(
								((pos[0] - min[0]) * factor).to_int_unchecked::<i32>(),
								((pos[1] - min[1]) * factor).to_int_unchecked::<i32>(),
							)
						}
					},
				],
				Into::<ShapeStyle>::into(&link_color).filled(),
			))
			.unwrap();
		}
	}

	let mut max_swinging = 0.0;
	let mut max_traction = 0.0;
	layout
		.speeds
		.iter()
		.zip(layout.old_speeds.iter())
		.for_each(|(speed, old_speed)| {
			let swinging = (speed[0] - old_speed[0]).powi(2) + (speed[1] - old_speed[1]).powi(2);
			let traction = (speed[0] + old_speed[0]).powi(2) + (speed[1] + old_speed[1]).powi(2);
			if swinging > max_swinging {
				max_swinging = swinging;
			}
			if traction > max_traction {
				max_traction = traction;
			}
		});

	for (pos, speed, old_speed) in izip!(
		layout.points.iter(),
		layout.speeds.iter(),
		layout.old_speeds.iter()
	) {
		let swinging = (speed[0] - old_speed[0]).powi(2) + (speed[1] - old_speed[1]).powi(2);
		let traction = (speed[0] + old_speed[0]).powi(2) + (speed[1] + old_speed[1]).powi(2);

		root.draw(&Circle::new(
			unsafe {
				(
					((pos[0] - min[0]) * factor).to_int_unchecked::<i32>(),
					((pos[1] - min[1]) * factor).to_int_unchecked::<i32>(),
				)
			},
			2,
			Into::<ShapeStyle>::into(&RGBColor(
				((swinging / max_swinging).powf(0.1) * 255.).min(255.) as u8,
				0,
				((traction / max_traction).powf(0.1) * 255.).min(255.) as u8,
			))
			.filled(),
		))
		.unwrap();
	}

	std::mem::drop(layout);

	println!("size:  {:?}", graph_size);
	println!("scale: {}", factor);
}
