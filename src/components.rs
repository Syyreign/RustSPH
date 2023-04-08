use bevy::prelude::*;

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

#[derive(Resource)]
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

#[derive(Resource)]
pub(crate) struct CollisionParameters {
    pub(crate) top_point: Vec3,
    pub(crate) bottom_point: Vec3,
}

impl Default for CollisionParameters{
    fn default() -> Self {
        CollisionParameters { 
            top_point: Vec3::ONE * 2.0,
            bottom_point: Vec3::NEG_ONE * 2.0,
        }
    }
}