use std::f32::consts::PI;

use bevy::app::App;
use bevy::prelude::*;
use bevy::time::FixedTimestep;

#[derive(Debug, Hash, PartialEq, Eq, Clone, StageLabel)]
struct FixedUpdateStage;

const DELTA_TIME: f64 = 0.01;
const PARTICLE_RADIUS: f32 = 0.75;
const SMOOTHING_RADIUS: f32 = 0.75;
const PRESSURE_CONSTANT: f32 = 20.0;
const REFERENCE_DENSITY: f32 = 20.0;

fn main(){
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugin(HelloPlugin)
        .run();
}

pub struct HelloPlugin;

impl Plugin for HelloPlugin {
    fn build(&self, app: &mut App) {
        app.add_startup_system(setup_scene)
            .add_startup_system(generate_fluid)
            .add_stage_after(
                CoreStage::Update,
                FixedUpdateStage,
                SystemStage::parallel()
                    .with_run_criteria(FixedTimestep::step(DELTA_TIME))
                    .with_system(clean_particle)
                    .with_system(update_pressure)
                    .with_system(update_acceleration)
                    .with_system(integrate),
            )
            .run();
    }
}

#[derive(Component)]
struct Particle;

#[derive(Component, Default)]
struct Mass(f32);
#[derive(Component, Default)]
struct Density(f32);
#[derive(Component, Default)]
struct Pressure(f32);
#[derive(Component, Default)]
struct Acceleration(Vec3);
#[derive(Component, Default)]
struct LastPos(Vec3);

#[derive(Bundle, Default)]
struct ParticleBundle{
    pbr: PbrBundle,
    mass: Mass,
    density: Density,
    pressure: Pressure,
    acceleration: Acceleration,
    last_pos: LastPos,
}

fn setup_scene(    
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
){
    commands.spawn(PointLightBundle {
        point_light: PointLight {
            intensity: 1500.0,
            shadows_enabled: true,
            ..default()
        },
        transform: Transform::from_xyz(4.0, 8.0, 4.0),
        ..default()
    });

    commands.spawn(Camera3dBundle {
        transform: Transform::from_xyz(-2.0, 2.5, 5.0).looking_at(Vec3::ZERO, Vec3::Y),
        ..default()
    });
}

fn generate_fluid(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let mesh = meshes.add(Mesh::from(shape::Icosphere {
        radius: 1.0,
        subdivisions: 3,
    }));

    for i in 0..4{

        let position = Vec3::new(
            0.0 + i as f32,
            0.5,
            0.0,
        );

        commands.spawn((ParticleBundle {
            pbr: PbrBundle {
                mesh: mesh.clone(),
                transform: Transform {
                    translation: position,
                    ..default()
                },
                material: materials.add(Color::rgb(0.8, 0.7, 0.6).into()),
                ..default()
            },
            mass: Mass(1.0),
            density: Density(0.0),
            pressure: Pressure(0.0),
            acceleration: Acceleration(Vec3{ x: 0.0, y: 0.0, z: 0.0}),
            last_pos: LastPos(position),
        },
        Particle,
        ));
    }

}

fn update_pressure(mut query: Query<(&Mass, &GlobalTransform, &mut Density, &mut Pressure), With<Particle>>) {
    let mut iter = query.iter_combinations_mut();

    // Set the density of each particle
    while let Some([(Mass(m1), transform1, mut density1, pressure1), 
        (Mass(m2), transform2, mut density2, pressure2)]) =
        iter.fetch_next()
    {
        let delta = transform2.translation() - transform1.translation();
        let distance_sq: f32 = delta.length_squared();
        if(distance_sq > PARTICLE_RADIUS){
            continue;
        }
        println!("trans1: {} trans 2: {} distance: {}", transform1.translation(), transform2.translation(), distance_sq);
        let poly6 = get_poly6_smoothing(distance_sq);
        density1.0 += m1 * poly6;
        density2.0 += m2 * poly6;
        println!("poly6: {} density1: {} density2: {}", poly6, density1.0, density2.0);
        let x = 5;

    }

    // Set the pressure of the particle
    for (_mass, _transform, density, mut pressure) in &mut query {
        pressure.0 = PRESSURE_CONSTANT * (density.0 - REFERENCE_DENSITY);
    }

    println!("");
}

fn update_acceleration(mut query: Query<(&Mass, &GlobalTransform, &Density, &Pressure, &mut Acceleration), With<Particle>>){
    let mut iter = query.iter_combinations_mut();

    // Set the density of each particle
    while let Some([(Mass(m1), transform1, density1, pressure1, mut accel1), 
        (Mass(m2), transform2, density2, pressure2, mut accel2)]) =
        iter.fetch_next()
    {
        let rij = transform1.translation() - transform2.translation();
        let rji = transform2.translation() - transform1.translation();
        let distance_sq: f32 = rij.length_squared();
        if distance_sq > PARTICLE_RADIUS {
            continue;
        }
        
        accel1.0 -= ((m2/m1) * ((pressure1.0 + pressure2.0) / (2.0 * density1.0 * density2.0))) 
            * get_spiky_smoothing(distance_sq) * rij;
        accel2.0 -= ((m1/m2) * ((pressure1.0 + pressure2.0) / (2.0 * density1.0 * density2.0)))
            * get_spiky_smoothing(distance_sq) * rji;
    }

    // Add gravity
    for (_mass, _transform, _density, _pressure, mut accel) in &mut query {
        accel.0 += Vec3{ x: 0.0, y: -9.8, z: 0.0};
    }
}

///
/// Input. distance_sq: The distance between two given particles
/// Ouput. The poly6 smoothing kernel
/// 
fn get_poly6_smoothing(distance_sq: f32) -> f32{
    let mut poly6 = ((SMOOTHING_RADIUS * SMOOTHING_RADIUS) - distance_sq);
    poly6 = f32::powf(poly6, 3.0);
    poly6 = (315.0 / (64.0 * PI * f32::powf(SMOOTHING_RADIUS, 9.0))) * poly6;

    poly6
}

///
/// Input. distance_sq: The distance between two given particles
/// Ouput. The poly6 smoothing kernel
/// 
fn get_spiky_smoothing(distance_sq: f32) -> f32{
    let mut spiky = f32::powf((SMOOTHING_RADIUS - f32::sqrt(distance_sq)), 2.0);
    spiky = f32::powf(spiky, 3.0);
    spiky = (-45.0 / (PI * f32::powf(SMOOTHING_RADIUS, 6.0))) * spiky;

    spiky
}

fn integrate(mut query: Query<(&mut Acceleration, &mut Transform, &mut LastPos)>) {
    let dt_sq = (DELTA_TIME * DELTA_TIME) as f32;
    for (mut acceleration, mut transform, mut last_pos) in &mut query {
        // verlet integration
        // x(t+dt) = 2x(t) - x(t-dt) + a(t)dt^2 + O(dt^4)

        let new_pos =
            transform.translation + transform.translation - last_pos.0 + acceleration.0 * dt_sq;
        acceleration.0 = Vec3::ZERO;
        last_pos.0 = transform.translation;
        transform.translation = new_pos;
    }
}

fn clean_particle(mut query: Query<(&mut Acceleration, &mut Density, &mut Pressure)>){
    for (mut acceleration, mut density, mut pressure) in &mut query {
        acceleration.0 = Vec3::ZERO;
        density.0 = 0.0;
        pressure.0 = 0.0;
    }
}