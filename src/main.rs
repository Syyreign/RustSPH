use std::f32::consts::PI;

use bevy::app::App;
use bevy::prelude::*;

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
struct FixedUpdateStage;

const DELTA_TIME: f32 = 0.001;
const SMOOTHING_RADIUS: f32 = 0.01;
const PRESSURE_CONSTANT: f32 = 20.0;
const REFERENCE_DENSITY: f32 = 20.0;
const MAX_ACCELERATION: f32 = 100.0;
const MAX_VELOCITY: f32 = 100.0;
const NEAR_ZERO: f32 = 0.0000001;

fn main(){
    App::new()
        .add_plugins(DefaultPlugins)
        .add_startup_system(setup_scene)
        .add_startup_system(generate_fluid)
        .insert_resource(FixedTime::new_from_secs(DELTA_TIME))
        .add_systems((
            clean_particle, 
            update_pressure, 
            //update_acceleration, 
            //plane_collision, 
            integrate,
            plane_collision
        ).in_schedule(CoreSchedule::FixedUpdate))
        .run();
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
struct Velocity(Vec3);
#[derive(Component, Default)]
struct LastPos(Vec3);

// Components for collision plane
#[derive(Component)]
struct Plane;
#[derive(Component, Default)]
struct Point(Vec3);
#[derive(Component, Default)]
struct Normal(Vec3);

#[derive(Bundle, Default)]
struct ParticleBundle{
    pbr: PbrBundle,
    mass: Mass,
    density: Density,
    pressure: Pressure,
    acceleration: Acceleration,
    velocity: Velocity,
    last_pos: LastPos,
}

#[derive(Bundle, Default)]
struct PlaneBundle{
    point: Point,
    normal: Normal,
}

fn setup_scene(    
    mut commands: Commands,
    mut _meshes: ResMut<Assets<Mesh>>,
    mut _materials: ResMut<Assets<StandardMaterial>>,
){
    commands.spawn(PointLightBundle {
        point_light: PointLight {
            intensity: 2500.0,
            shadows_enabled: true,
            ..default()
        },
        transform: Transform::from_xyz(4.0, 8.0, 4.0),
        ..default()
    });

    let points = 
        [Vec3{ x: 0.0, y: -80.0, z: 0.0}, 
        Vec3{ x: 0.0, y: 80.0, z: 0.0}, 
        Vec3{ x: -20.0, y: 0.0, z: 0.0}, 
        Vec3{ x: 20.0, y: 0.0, z: 0.0}, 
        Vec3{ x: 0.0, y: 0.0, z: -20.0}, 
        Vec3{ x: 0.0, y: 0.0, z: 20.0}];

    let normals = [Vec3{ x: 0.0, y: 1.0, z: 0.0}, 
        Vec3{ x: 0.0, y: -1.0, z: 0.0}, 
        Vec3{ x: 1.0, y: 0.0, z: 0.0}, 
        Vec3{ x: -1.0, y: 0.0, z: 0.0}, 
        Vec3{ x: 0.0, y: 0.0, z: 1.0}, 
        Vec3{ x: 0.0, y: 0.0, z: -1.0}];

    for i in 0..points.len(){
        commands.spawn((PlaneBundle {
            point: Point(points[i]),
            normal: Normal(normals[i]),
        },
        Plane
        ));
    }

    commands.spawn(Camera3dBundle {
        transform: Transform::from_xyz(150.0, 70.0, -150.0).looking_at(Vec3::ZERO, Vec3::Y),
        ..default()
    });

    
}

fn generate_fluid(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {

    let mesh = meshes.add(
        Mesh::try_from(shape::Icosphere {
            radius: 0.5,
            subdivisions: 3,
        })
        .unwrap(),
    );

    for i in 0..2{
        for j in 0..2{
            for k in 0..2{

                let position = Vec3::new(
                    0.0 + (i as f32 * 0.5) + (j as f32 * 0.6),
                    0.0 + (j as f32 * 2.0),
                    0.0 + (k as f32 * 0.5) + (j as f32 * 0.6),
                );

                commands.spawn((ParticleBundle {
                    pbr: PbrBundle {
                        transform: Transform {
                            translation: position,
                            ..default()
                        },
                        mesh: mesh.clone(),
                        material: materials.add(Color::rgb(0.1, 0.1, 0.8).into()),
                        ..default()
                    },
                    mass: Mass(10.0),
                    density: Density(0.0),
                    pressure: Pressure(0.0),
                    acceleration: Acceleration(Vec3::ZERO),
                    velocity: Velocity(Vec3::ZERO),
                    last_pos: LastPos(position),
                },
                Particle,
                ));
            }
        }
    }
    

}

fn update_pressure(mut query: Query<(&Mass, &GlobalTransform, &mut Density, &mut Pressure, &mut Acceleration), With<Particle>>) {
    
    {
        let mut iter = query.iter_combinations_mut();

        // Set the density of each particle
        while let Some([(Mass(m1), transform1, mut density1, _pressure1, mut _accel1), 
            (Mass(m2), transform2, mut density2, _pressure2, mut _accel2)]) =
            iter.fetch_next()
        {
            let delta = transform2.translation() - transform1.translation();
            let distance_sq: f32 = delta.length_squared();

            let poly6 = get_poly6_smoothing(distance_sq);
            density1.0 += *m1 * poly6;
            density2.0 += *m2 * poly6;
        }
    }

    // Set the pressure of the particle
    for (_mass, _transform, mut density, mut pressure, mut _accel) in &mut query {
        if density.0 < REFERENCE_DENSITY{
            density.0 = REFERENCE_DENSITY;
        }
        
        pressure.0 = PRESSURE_CONSTANT * (density.0 - REFERENCE_DENSITY);
        let x =5;
    }

    {
        let mut iter = query.iter_combinations_mut();

        // Set the density of each particle
        while let Some([(Mass(m1), transform1, density1, pressure1, mut accel1), 
            (Mass(m2), transform2, density2, pressure2, mut accel2)]) =
            iter.fetch_next()
        {
            let rij = transform1.translation() - transform2.translation();
            let rji = transform2.translation() - transform1.translation();
            let distance_sq: f32 = rij.length_squared();

            let spiky = get_spiky_smoothing(distance_sq);

            // Acceleration
            if density1.0.abs() > NEAR_ZERO && density2.0.abs() > NEAR_ZERO {
                accel1.0 += -(( *m2/ *m1) * ((pressure1.0 + pressure2.0) / (2.0 * density1.0 * density2.0))) 
                    * spiky * rij.normalize();
                accel2.0 += -(( *m1/ *m2) * ((pressure1.0 + pressure2.0) / (2.0 * density1.0 * density2.0)))
                * spiky * rji.normalize();
            }

            // Viscocity Acceleration
            let epsillon = 0.018;

            if density2.0 > NEAR_ZERO { 
                let av1 = epsillon * (m2/m1) * (1.0 / density2.0) * (rji.normalize()) * get_viscosity_smoothing(distance_sq);
                accel1.0 += av1;

            }
            if density1.0 > NEAR_ZERO {
                let av2 = epsillon * (m1/m2) * (1.0 / density1.0) * (rij.normalize()) * get_viscosity_smoothing(distance_sq);
                accel2.0 += av2;
            }
        }
    }

    // Add gravity and clamp
    for (_mass, _transform, _density, _pressure, mut accel) in &mut query {
        accel.0 += Vec3{ x: 0.0, y: -9.8, z: 0.0};

        if accel.0.length() > MAX_ACCELERATION {
            accel.0 = accel.0.normalize() * MAX_ACCELERATION;
        }
    }

    //println!("Sim");
}

fn accel(mut query: Query<(&Mass, &GlobalTransform, &Density, &Pressure, &mut Acceleration), With<Particle>>){
    for (_mass, _transform, density, _pressure, mut accel) in &mut query {
        println!("dens: {}", density.0);
    }
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

        let spiky =get_spiky_smoothing(distance_sq);

        // Acceleration
        if density1.0.abs() > NEAR_ZERO && density2.0.abs() > NEAR_ZERO {
            accel1.0 += (( *m2/ *m1) * ((pressure1.0 + pressure2.0) / (2.0 * density1.0 * density2.0))) 
                * spiky * rij.normalize();
            accel2.0 += (( *m1/ *m2) * ((pressure1.0 + pressure2.0) / (2.0 * density1.0 * density2.0)))
            * spiky * rji.normalize();
        }

        // Viscocity Acceleration
        let epsillon = 0.018;

        if density2.0 > NEAR_ZERO { 
            let av1 = epsillon * (m2/m1) * (1.0 / density2.0) * (rji.normalize()) * get_viscosity_smoothing(distance_sq);
            accel1.0 += av1;

        }
        if density1.0 > NEAR_ZERO {
            let av2 = epsillon * (m1/m2) * (1.0 / density1.0) * (rij.normalize()) * get_viscosity_smoothing(distance_sq);
            accel2.0 += av2;
        }
    }

    // Add gravity and clamp
    for (_mass, _transform, _density, _pressure, mut accel) in &mut query {
        accel.0 += Vec3{ x: 0.0, y: -9.8, z: 0.0};

        if accel.0.length() > MAX_ACCELERATION {
            accel.0 = accel.0.normalize() * MAX_ACCELERATION;
        }
    }


}

fn plane_collision(
    mut particles: Query<(&mut Acceleration, &mut Velocity, &mut Transform), With<Particle>>,
    mut planes: Query<(&Point, &Normal), With<Plane>>,
){
    for (mut _acceleration, mut velocity, transform) in &mut particles {
        for (point, normal) in &mut planes{

            // Point is colliding if < 0.0
            if (transform.translation - point.0).dot(normal.0) < 0.0 && velocity.0.dot(normal.0) < 0.0{
                velocity.0 = -((0.75 * normal.0).dot(velocity.0)) * (normal.0);
                //velocity.0 += normal.0 * 0.5;
            }
        }
        
    }
}

///
/// Input. distance_sq: The distance between two given particles
/// Ouput. The poly6 smoothing kernel
/// 
fn get_poly6_smoothing(distance_sq: f32) -> f32{
    let mut poly6 = (SMOOTHING_RADIUS * SMOOTHING_RADIUS) - distance_sq;

    poly6 = f32::powf(poly6, 3.0);
    poly6 = (315.0 / (64.0 * PI * f32::powf(SMOOTHING_RADIUS, 9.0))) * poly6;

    poly6
}

///
/// Input. distance_sq: The distance between two given particles
/// Ouput. The poly6 smoothing kernel
/// 
fn get_spiky_smoothing(distance_sq: f32) -> f32{
    let mut spiky = f32::powf(SMOOTHING_RADIUS - f32::sqrt(distance_sq), 2.0);
    spiky = (-45.0 / (PI * f32::powf(SMOOTHING_RADIUS, 6.0))) * spiky;

    spiky
}

fn get_viscosity_smoothing(distance_sq: f32) -> f32{
    (45.0 / (PI * f32::powf(SMOOTHING_RADIUS, 6.0))) * (SMOOTHING_RADIUS - f32::sqrt(distance_sq))
}

fn integrate(mut query: Query<(&mut Acceleration, &mut Velocity, &mut Transform, &mut LastPos)>) {
    // let dt_sq = (DELTA_TIME * DELTA_TIME) as f32;
    // for (mut acceleration, mut velocity, mut transform, mut last_pos) in &mut query {
    //     // verlet integration
    //     // x(t+dt) = 2x(t) - x(t-dt) + a(t)dt^2 + O(dt^4)
    //     println!("Integrate. accel: {}", acceleration.0);
    //     //velocity.0 = veloctiy.0 + 

    //     let new_pos =
    //         (transform.translation + transform.translation) - last_pos.0 + (acceleration.0 * dt_sq);
    //     acceleration.0 = Vec3::ZERO;
    //     last_pos.0 = transform.translation;
    //     transform.translation = new_pos;
    // }

    for (mut acceleration, mut velocity, mut transform, mut last_pos) in &mut query {
        velocity.0 = velocity.0 + (acceleration.0 * DELTA_TIME as f32);
        if velocity.0.length() > MAX_VELOCITY {
            velocity.0 = velocity.0.normalize() * MAX_VELOCITY;
        }

        transform.translation = transform.translation + (velocity.0 * DELTA_TIME as f32);
        acceleration.0 = Vec3::ZERO;
    }
    //println!("Integrate");
}

fn clean_particle(mut query: Query<(&mut Acceleration, &mut Velocity, &mut Density, &mut Pressure)>){
    for (mut _acceleration, mut _velocity, mut density, mut pressure) in &mut query {
        //acceleration.0 = Vec3::ZERO;
        //velocity.0 = Vec3::ZERO;
        density.0 = 0.0;
        pressure.0 = 0.0;
    }
    let x = 5;
}