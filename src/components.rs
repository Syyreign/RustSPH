use bevy::prelude::*;

#[derive(Component)]
pub(crate) struct GameCamera;

#[derive(Component)]
pub(crate) struct Particle;
#[derive(Component, Default)]
pub(crate) struct Mass(pub(crate) f32);
#[derive(Component, Default)]
pub(crate) struct Density(pub(crate) f32);
#[derive(Component, Default)]
pub(crate) struct Pressure(pub(crate) f32);
#[derive(Component, Default)]
pub(crate) struct Acceleration(pub(crate) Vec3);
#[derive(Component, Default)]
pub(crate) struct Velocity(pub(crate) Vec3);
#[derive(Component, Default)]
pub(crate) struct LastPos(pub(crate) Vec3);

// Components for collision plane
#[derive(Component)]
pub(crate) struct Plane;
#[derive(Component, Default)]
pub(crate) struct Point(pub(crate) Vec3);
#[derive(Component, Default)]
pub(crate) struct Normal(pub(crate) Vec3);

#[derive(Bundle, Default)]
pub(crate) struct ParticleBundle{
    pub(crate) pbr: PbrBundle,
    pub(crate) mass: Mass,
    pub(crate) density: Density,
    pub(crate) pressure: Pressure,
    pub(crate) acceleration: Acceleration,
    pub(crate) velocity: Velocity,
    pub(crate) last_pos: LastPos,
}

#[derive(Bundle, Default)]
pub(crate) struct PlaneBundle{
    pub(crate) point: Point,
    pub(crate) normal: Normal,
}

#[derive(Resource, Clone)]
pub(crate) struct FluidParameters {
    pub(crate) delta_time: f32,
    pub(crate) smoothing_radius: f32,
    pub(crate) pressure_constant: f32,
    pub(crate) reference_density: f32,
    pub(crate) max_acceleration: f32,
    pub(crate) max_velocity: f32,
    pub(crate) mass: f32,
    pub(crate) viscosity_coef: f32,
    pub(crate) gravity: f32,
    pub(crate) near_zero: f32,
}

impl Default for FluidParameters{
    fn default() -> Self {
        FluidParameters { 
            delta_time: 0.01,
            smoothing_radius: 0.5,
            pressure_constant: 25.0,
            reference_density: 15.0,
            max_acceleration: 100.0,
            max_velocity: 100.0,
            mass: 1.0,
            viscosity_coef: 0.02, 
            gravity: -9.8,
            near_zero: 0.0000001,
        }
    }
}

#[derive(Resource, Clone)]
pub(crate) struct SimulationParameters {
    pub(crate) top_point: Vec3,
    pub(crate) bottom_point: Vec3,
    pub(crate) num_particles: Vec3,
    pub(crate) particle_distance: f32,
    pub(crate) particle_origin: Vec3,
    pub(crate) camera_position: Vec3,
    pub(crate) camera_target: Vec3,
    pub(crate) planet_gravity: bool,
}

impl Default for SimulationParameters{
    fn default() -> Self {
        SimulationParameters { 
            top_point: Vec3::new(2.0, 2.0, 2.0),
            bottom_point: Vec3::new(-2.0, -2.0 ,-2.0),
            num_particles: Vec3::new(6.0, 6.0 ,6.0),
            particle_distance: 0.5,
            particle_origin: Vec3::new(0.0, 0.0 ,0.0),
            camera_position: Vec3::new(6.0, 3.0, 6.0),
            camera_target: Vec3::new(0.0, 0.0, 0.0),
            planet_gravity: false,
        }
    }
}

#[derive(Clone)]
pub(crate) struct Demo{
    pub(crate) fluid_parameters: FluidParameters,
    pub(crate) simulation_parameters: SimulationParameters,
}

impl Default for Demo{
    fn default() -> Self {
        Demo{
            fluid_parameters: FluidParameters::default(),
            simulation_parameters: SimulationParameters::default(),
        }
    }
}

pub(crate) const DEMOS: [Demo;4] = 
[
    Demo{
        fluid_parameters: FluidParameters { 
            delta_time: 0.01,
            smoothing_radius: 0.5,
            pressure_constant: 25.0,
            reference_density: 15.0,
            max_acceleration: 100.0,
            max_velocity: 100.0,
            mass: 1.0,
            viscosity_coef: 0.02, 
            gravity: -9.8,
            near_zero: 0.0000001,
        },
        simulation_parameters: SimulationParameters { 
            top_point: Vec3::new(2.0, 2.0, 2.0),
            bottom_point: Vec3::new(-2.0, -2.0 ,-2.0),
            num_particles: Vec3::new(6.0, 6.0 ,6.0),
            particle_distance: 0.5,
            particle_origin: Vec3::new(0.0, 0.0 ,0.0),
            camera_position: Vec3::new(6.0, 3.0, 6.0),
            camera_target: Vec3::ZERO,
            planet_gravity: false,
        }
    },
    Demo{
        fluid_parameters: FluidParameters { 
            delta_time: 0.015,
            smoothing_radius: 0.5,
            pressure_constant: 25.0,
            reference_density: 15.0,
            max_acceleration: 150.0,
            max_velocity: 150.0,
            mass: 1.0,
            viscosity_coef: 0.1, 
            gravity: -9.8,
            near_zero: 0.0000001,
        },
        simulation_parameters: SimulationParameters { 
            top_point: Vec3::new(2.0, 2.0, 2.0),
            bottom_point: Vec3::new(-2.0, -2.0 ,-2.0),
            num_particles: Vec3::new(7.0, 9.0 ,7.0),
            particle_distance: 0.5,
            particle_origin: Vec3::new(0.0, 0.0 ,0.0),
            camera_position: Vec3::new(6.0, 3.0, 6.0),
            camera_target: Vec3::ZERO,
            planet_gravity: false,
        }
    },
    Demo{
        fluid_parameters: FluidParameters { 
            delta_time: 0.015,
            smoothing_radius: 0.5,
            pressure_constant: 25.0,
            reference_density: 15.0,
            max_acceleration: 150.0,
            max_velocity: 150.0,
            mass: 1.0,
            viscosity_coef: 0.1, 
            gravity: -9.8,
            near_zero: 0.0000001,
        },
        simulation_parameters: SimulationParameters { 
            top_point: Vec3::new(3.0, 3.0, 3.0),
            bottom_point: Vec3::new(-3.0, -3.0 ,-3.0),
            num_particles: Vec3::new(7.0, 9.0 ,7.0),
            particle_distance: 0.5,
            particle_origin: Vec3::new(0.0, 0.0 ,0.0),
            camera_position: Vec3::new(6.0, 3.0, 6.0),
            camera_target: Vec3::ZERO,
            planet_gravity: true,
        }
    },
    Demo{
        fluid_parameters: FluidParameters { 
            delta_time: 0.015,
            smoothing_radius: 0.5,
            pressure_constant: 25.0,
            reference_density: 15.0,
            max_acceleration: 150.0,
            max_velocity: 150.0,
            mass: 1.0,
            viscosity_coef: 0.1, 
            gravity: -9.8,
            near_zero: 0.0000001,
        },
        simulation_parameters: SimulationParameters { 
            top_point: Vec3::new(0.5, 2.0, 2.0),
            bottom_point: Vec3::new(-0.5, -2.0 ,-2.0),
            num_particles: Vec3::new(6.0, 9.0 ,9.0),
            particle_distance: 0.3,
            particle_origin: Vec3::new(0.0, 0.0 ,0.0),
            camera_position: Vec3::new(16.0, 0.0, 7.0),
            camera_target: Vec3::new(0.0, 0.0, 7.0),
            planet_gravity: false,
        }
    }
];