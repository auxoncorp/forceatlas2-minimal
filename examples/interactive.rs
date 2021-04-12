use forceatlas2::*;
use itertools::izip;
use plotters::prelude::*;
use sdl2;
use sdl2::{event::Event, keyboard::Keycode, pixels::PixelFormatEnum};
use std::{
	io::BufRead,
	sync::{Arc, RwLock},
	thread,
	time::Duration,
};

const DRAW_LINKS: bool = true;
const SIZE: (u32, u32) = (800, 800);
const FRAMEDUR: u64 = 30;

pub fn run(image: Arc<RwLock<Vec<u8>>>) {
	let sdl_context = sdl2::init().unwrap();
	let video_subsystem = sdl_context.video().unwrap();

	let window = video_subsystem
		.window("ForceAtlas2", 800, 800)
		.resizable()
		.position_centered()
		.opengl()
		.build()
		.unwrap();

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
				_ => {}
			}
		}

		let mut texture = texture_creator
			.create_texture_streaming(PixelFormatEnum::RGB24, SIZE.0, SIZE.1)
			.unwrap();
		texture
			.update(None, &image.read().unwrap(), (SIZE.0 * 3u32) as usize)
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
		dimensions: 2,
		dissuade_hubs: false,
		ka: 0.01,
		kg: 0.001,
		kr: 0.002,
		lin_log: false,
		prevent_overlapping: None,
		strong_gravity: false,
	};

	let layout = Arc::new(RwLock::new(Layout::<f64>::from_graph(
		edges.clone(),
		Nodes::Degree(nodes),
		settings.clone(),
	)));
	let image = Arc::new(RwLock::new(vec![
		0u8;
		SIZE.0 as usize * SIZE.1 as usize * 3
	]));
	let computing = Arc::new(RwLock::new(false));
	let sleep = Arc::new(RwLock::new(50u64));
	{
		let layout = layout.write().unwrap();
		draw_graph(&layout, image.clone());
	}

	thread::spawn({
		let image = image.clone();
		let computing = computing.clone();
		let layout = layout.clone();
		let iters = iters.clone();
		move || loop {
			if *computing.read().unwrap() {
				let mut layout = layout.write().unwrap();
				draw_graph(&layout, image.clone());
				layout.iteration();
				iters.write().unwrap().add_assign(1);
			} else {
				thread::sleep(std::time::Duration::from_millis(*sleep.read().unwrap()));
			}
		}
	});

	thread::spawn(move || run(image));

	thread::spawn(move || loop {
		let cmd = scanrs::scanln();
		let mut args = cmd.split(' ');
		match args.next() {
			Some("q") => break,
			Some("h") => {
				println!("h         help");
				println!("q         quit");
				println!("r         reset");
				println!("s         start/stop");
				println!("ka <f64>  attraction");
				println!("kg <f64>  gravity");
				println!("kr <f64>  repulsion");
				println!("jt <f64>  jitter tolerance");
			}
			Some("s") => {
				let mut computing = computing.write().unwrap();
				*computing = !*computing;
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
				*layout = Layout::from_graph(edges.clone(), Nodes::Degree(nodes), settings.clone());
			}
			Some("p") => println!("ka={}  kg={}  kr={}", settings.ka, settings.kg, settings.kr,),
			Some("i") => println!("{}", iters.read().unwrap()),
			_ => println!("Unknown command"),
		}
	})
	.join()
	.unwrap();
}

fn draw_graph(layout: &Layout<f64>, image: Arc<RwLock<Vec<u8>>>) {
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
		let factors = (SIZE.0 as f64 / graph_size.0, SIZE.1 as f64 / graph_size.1);
		if factors.0 > factors.1 {
			min[0] -= (SIZE.0 as f64 / factors.1 - graph_size.0) / 2.0;
			factors.1
		} else {
			min[1] -= (SIZE.1 as f64 / factors.0 - graph_size.1) / 2.0;
			factors.0
		}
	};
	println!("size:  {:?}", graph_size);
	println!("scale: {}", factor);

	let mut image = image.write().unwrap();
	let root = BitMapBackend::with_buffer(&mut image, SIZE).into_drawing_area();
	root.fill(&WHITE).unwrap();

	if DRAW_LINKS {
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
				Into::<ShapeStyle>::into(&RGBColor(5, 5, 5).mix(0.05)).filled(),
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
}
