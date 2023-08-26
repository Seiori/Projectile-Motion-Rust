use nalgebra::{Vector3};
use std::sync::{Arc, Mutex};
use bmp::{Image, Pixel};
use scoped_threadpool::Pool;
use rand::Rng;

const WIDTH : usize = 1000;
const HEIGHT : usize = 1000;
const NUM_OF_PARTICLES : usize = 300000;
const CONE_ANGLE : f32 = 35.0;
const BLEND_COEFFICIENT : f32 = 0.1;

struct Paper {
    pixels: Arc<Mutex<Vec<Vec<Particle>>>>,
}

impl Paper {
    fn new() -> Self {
        let paper = Self {
            pixels: Arc::new(Mutex::new(vec![vec![Particle::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0)); WIDTH]; HEIGHT])),
        };

        // Change position value for all elements
        for x in 0..WIDTH {
            for y in 0..HEIGHT {
                let position_x = x as f32 * 0.001;
                let position_y = 0.0;
                let position_z = y as f32 * 0.001;

                paper.pixels.lock().unwrap()[x][y].position = Vector3::new(position_x, position_y, position_z);
            }
        }
        paper
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Particle {
    position: Vector3<f32>,
    colour: Vector3<f32>,
}

impl Particle {
    fn new(pos: Vector3<f32>, col: Vector3<f32>) -> Self {
        Self { 
            position: pos,
            colour: col,
        }
    }
}

struct ParticleSystem {
    particles: Vec<Particle>,
}

impl ParticleSystem {
    fn new(pos: Vector3<f32>, col: Vector3<f32>) -> Self {
        Self {
            particles: (0..NUM_OF_PARTICLES).map(|_| Particle::new(pos, col)).collect(),
        }
    }

    fn move_particles(particle: &mut Particle) {
        let mut rng = rand::thread_rng();
        let spray_cone = Vector3::new(
            rng.gen_range(-CONE_ANGLE..CONE_ANGLE),
            rng.gen_range(-CONE_ANGLE..CONE_ANGLE),
            rng.gen_range(-CONE_ANGLE..CONE_ANGLE),
        );  
        let spray_direction = (Vector3::new((WIDTH / 2) as f32, 0.0, (HEIGHT / 2) as f32) - particle.position).normalize();
        let rotated_direction = rotate_vector(spray_direction, spray_cone);
      
        let time_step = 0.01;
        let initial_velocity = 100.0;
        let gravity = Vector3::new(0.0, -9.8, 0.0); // Assuming gravity points downwards
        let drag_coefficient = 0.025; // Adjust this value based on the desired drag effect
      
        let mut velocity = rotated_direction * initial_velocity;

        // Apply drag
        let drag = -drag_coefficient * velocity;
        // Calculate net acceleration
        let net_acceleration = gravity + drag;

        while particle.position.y > 0.0 {
            // Update velocity using Euler's method
            velocity += net_acceleration * time_step;
            // Update position
            particle.position += velocity * time_step;
        }
    }
    
    fn run_simulation(&mut self) {
        // Create a Local Mutable Reference to the Particle System's Particles
        let particles = &mut self.particles;

        // Create a new Thread Pool
        let mut pool = Pool::new(6);

        pool.scoped(|scope| {
            for slice in particles.chunks_mut(NUM_OF_PARTICLES / 4) {
                scope.execute(move || thread_main_move(slice));
            }
        });
    }
}

fn main() {
    // Time the Simulation
    let start_time = std::time::SystemTime::now();

    // Generate a new Paper
    let paper = Paper::new();

    // Print the Time Taken to Move Particles
    println!("Time Taken to Initialize Paper Array Elements: {}ms", start_time.elapsed().unwrap().as_millis());

    // Generate Multiple Particle Systems
    let mut spray_cans = vec![
        ParticleSystem::new(Vector3::new((WIDTH + 50) as f32, 20.0, (HEIGHT / 5) as f32) ,Vector3::new(0.0, 0.0, 1.0)),
        ParticleSystem::new(Vector3::new(-50.0, 20.0, (HEIGHT / 5) as f32) ,Vector3::new(1.0, 0.0, 0.0)),
        ParticleSystem::new(Vector3::new((WIDTH / 2) as f32, 20.0, (HEIGHT + 50) as f32) ,Vector3::new(0.0, 1.0, 0.0)),
    ];

    // Print the Time Taken to Move Particles
    println!("Time Taken to Initialize All Spray Can Elements: {}ms", start_time.elapsed().unwrap().as_millis());

    // Run the Simulation for each Particle System
    let mut pool = Pool::new(3);
    pool.scoped(|scope| {
        for slice in spray_cans.chunks_mut(1) {
            scope.execute(move || thread_main_cans(slice));
        }
    });

    // Print the Time Taken to Move Particles
    println!("Time Taken to Calculate Particle Movements: {}ms", start_time.elapsed().unwrap().as_millis());
    
    let mut total_collided_particles = 0;

    // For each Particle System, check if the Particles have collided with the Paper's Pixels
    for spray_can in &spray_cans {
        for particle in &spray_can.particles {
            let x = particle.position.x.round() as usize;
            let y = particle.position.z.round() as usize;
    
            if x > 0 && x < WIDTH && y > 0 && y < HEIGHT {
                total_collided_particles += 1;

                // Find the Closest Pixel to the Particle
                let closest_pixel = paper.pixels.lock().unwrap()[x][y];
    
                // Change the colour of the Pixel
                paper.pixels.lock().unwrap()[x][y].colour = (1.0 - BLEND_COEFFICIENT) * closest_pixel.colour + BLEND_COEFFICIENT * particle.colour;
            }
        }
    }

    // Print the Time Taken to Check for Collisions Between the Paper and Particles, and to Blend the Colours
    println!("Time Taken to Check for Collisions Between the Paper and Particles, and to Blend the Colours: {}ms", start_time.elapsed().unwrap().as_millis());

    // Display the 2D Array of Paper as an Image
    let img = {
        
        let mut img = Image::new(WIDTH as u32, HEIGHT as u32);
        for x in 0..WIDTH {
            for y in 0..HEIGHT {
                let pixel = paper.pixels.lock().unwrap()[x][y];
                img.set_pixel(x as u32, y as u32, Pixel::new(
                    (pixel.colour.x * 255.0) as u8,
                    (pixel.colour.y * 255.0) as u8,
                    (pixel.colour.z * 255.0) as u8,
                ));
            }
        }
        img
    };
    // Save the Image
    img.save("output.bmp").unwrap();

    // Print the Total Time Taken
    println!("Time Taken to Save the Produced Image: {}ms", start_time.elapsed().unwrap().as_millis());

    // Print Simulation Statistics
    println!("Number of Particles Created: {}", NUM_OF_PARTICLES * 3);
    println!("Number of Particles that Collided with the Paper: {}", total_collided_particles);
    println!("Number of Particles that didn't Collide with the Paper: {}", NUM_OF_PARTICLES * 3 - total_collided_particles);
}

fn thread_main_move(list: &mut [Particle]) {
    for i in 0..list.len() {
        ParticleSystem::move_particles(&mut list[i]);
    }
}

fn thread_main_cans(list: &mut [ParticleSystem]) {
    for i in 0..list.len() {
        list[i].run_simulation();
    }
}

fn rotate_vector(vector: Vector3<f32>, angles: Vector3<f32>) -> Vector3<f32> {
    let angle_x = angles.x.to_radians();
    let angle_y = angles.y.to_radians();
    let angle_z = angles.z.to_radians();

    let rotation_matrix = nalgebra::Rotation3::from_euler_angles(angle_x, angle_y, angle_z);

    rotation_matrix * vector
}