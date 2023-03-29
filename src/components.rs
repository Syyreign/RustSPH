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