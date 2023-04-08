use std::f32::consts::PI;

use bevy::app::App;
use bevy::prelude::*;
use bevy_egui::{egui, EguiContexts, EguiPlugin};
use bevy_prototype_debug_lines::*;

use bevy::diagnostic::{FrameTimeDiagnosticsPlugin, LogDiagnosticsPlugin};

use components::*;
use rand::Rng;

mod components;

#[derive(Debug, Hash, PartialEq, Eq, Clone)]
struct FixedUpdateStage;

fn main(){
    App::new()
        .add_plugins(DefaultPlugins)
        .add_plugin(EguiPlugin)
        .add_plugin(DebugLinesPlugin::default())
        .add_startup_system(setup_scene)
        .add_startup_system(generate_fluid)
        .add_plugin(LogDiagnosticsPlugin::default())
        .add_plugin(FrameTimeDiagnosticsPlugin::default())
        .insert_resource(ClearColor(Color::rgb(1.0, 1.0, 1.0)))
        .init_resource::<FluidParameters>()
        .init_resource::<CollisionParameters>()
        .add_system(spawn_particle,)
        .add_system(ui_system)
        .add_systems((
            draw_bounds,
            clean_particle, 
            update_pressure
                .after(clean_particle), 
            update_acceleration
                .after(update_pressure), 
            plane_repulsion
                .after(update_acceleration), 
            integrate
                .after(update_pressure),
            plane_collision
                .after(integrate)
        ).in_schedule(CoreSchedule::FixedUpdate))
        .run();
}

fn ui_system(mut contexts: EguiContexts, mut fluid_params: ResMut<FluidParameters>, mut coll_params: ResMut<CollisionParameters>) {
    egui::Window::new("Controls").show(contexts.ctx_mut(), |ui| {
        ui.label("Fluid Controls");

        ui.add(
            egui::Slider::new(&mut fluid_params.delta_time, 0.0..=0.5)
            .text("Delta Time")
            .logarithmic(true),
        );
        ui.add(
            egui::Slider::new(&mut fluid_params.smoothing_radius, 0.0..=10.)
            .text("Smoothing Radius"),
        );
        ui.add(
            egui::Slider::new(&mut fluid_params.pressure_constant, 0.0..=100.)
            .text("Pressure Constant"),
        );
        ui.add(
            egui::Slider::new(&mut fluid_params.reference_density, 0.0..=100.)
            .text("Reference Density"),
        );
        ui.add(
            egui::Slider::new(&mut fluid_params.max_acceleration, 0.0..=1000.)
            .text("Max Acceleration"),
        );
        ui.add(
            egui::Slider::new(&mut fluid_params.max_velocity, 0.0..=1000.)
            .text("Max Velocity"),
        );
        ui.add(
            egui::Slider::new(&mut fluid_params.mass, 0.0..=100.)
            .text("Mass"),
        );
        ui.add(
            egui::Slider::new(&mut fluid_params.viscosity_coef, 0.0..=10.)
            .text("Viscosity Coef"),
        );

        // Collision. TODO make this not repeat
        egui::Grid::new("top_point")
            .num_columns(3)
            .spacing([4.0, 4.0])
            .striped(true)
            .show(ui, |ui| 
        {
            ui.add(
                egui::DragValue::new(&mut coll_params.top_point[0])
                .clamp_range(0.1..=10.)
                .speed(0.1),
            );
            ui.add(
                egui::DragValue::new(&mut coll_params.top_point[1])
                .clamp_range(0.1..=10.)
                .speed(0.1),
            );
            ui.add(
                egui::DragValue::new(&mut coll_params.top_point[2])
                .clamp_range(0.1..=10.)
                .speed(0.1),
            );
            ui.add(
                egui::Label::new("Top Point")
            );
            ui.end_row();
        });
        
        egui::Grid::new("bottom_point")
            .num_columns(3)
            .spacing([4.0, 4.0])
            .striped(true)
            .show(ui, |ui| 
        {
            ui.add(
                egui::DragValue::new(&mut coll_params.bottom_point[0])
                .clamp_range(-0.1..=-10.)
                .speed(0.1),
            );
            ui.add(
                egui::DragValue::new(&mut coll_params.bottom_point[1])
                .clamp_range(-0.1..=-10.)
                .speed(0.1),
            );
            ui.add(
                egui::DragValue::new(&mut coll_params.bottom_point[2])
                .clamp_range(-0.1..=-10.)
                .speed(0.1),
            );
            ui.add(
                egui::Label::new("Bottom Point")
            );
            ui.end_row();
        });
    });
}

fn draw_bounds(
    coll_params: Res<CollisionParameters>,
    mut lines: ResMut<DebugLines>,
    mut query: Query<&mut components::Point, With<components::Plane>>
){

    let mut iter = query.iter_combinations();

    while let Some([point1, 
        point2,
        point3]) =
        iter.fetch_next()
    {
        let corner = point1.0 + point2.0 + point3.0;

        //If an axis is zero continue
        if (corner[0] * corner[1] * corner[2]).abs() < 0.00001 {
            continue;
        }

        for i in 0..3{
            let mut mult = Vec3::ONE;
            mult[i] = mult[i] * 0.9; 

            lines.line_colored(corner, corner * mult, 0.0, Color::rgb(0.0, 0.0, 0.0));
        }
    }

    // Gross code to update the line position of the bounds
    for (mut point) in &mut query{
        if(point.0[0] != 0.0){
            if(point.0[0] < 0.0){
                point.0[0] = coll_params.bottom_point[0];
                continue;
            }
            point.0[0] = coll_params.top_point[0];
        }

        if(point.0[1] != 0.0){
            if(point.0[1] < 0.0){
                point.0[1] = coll_params.bottom_point[1];
                continue;
            }
            point.0[1] = coll_params.top_point[1];
        }

        if(point.0[2] != 0.0){
            if(point.0[2] < 0.0){
                point.0[2] = coll_params.bottom_point[2];
                continue;
            }
            point.0[2] = coll_params.top_point[2];
        }
    }

}

fn spawn_particle(
    params: Res<FluidParameters>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    keyboard_input: Res<Input<KeyCode>>,
){
    if keyboard_input.just_pressed(KeyCode::P) {

        for _ in 0..10{
            let mesh = meshes.add(
                Mesh::try_from(shape::Icosphere {
                    radius: 0.1,
                    subdivisions: 3,
                })
                .unwrap(),
            );

            let mut rng = rand::thread_rng();
            let position = Vec3::new(rng.gen_range(-0.5..0.5), 1.0, rng.gen_range(-0.5..0.5));

            commands.spawn((components::ParticleBundle {
                pbr: PbrBundle {
                    transform: Transform {
                        translation: position,
                        ..default()
                    },
                    mesh: mesh.clone(),
                    material: materials.add(Color::rgb(0.4, 0.98, 1.0).into()),
                    ..default()
                },
                mass: components::Mass(params.mass),
                density: components::Density(0.0),
                pressure: components::Pressure(0.0),
                acceleration: components::Acceleration(Vec3::ZERO),
                velocity: components::Velocity(Vec3::ZERO),
                last_pos: components::LastPos(position),
            },
            components::Particle,
            ));
        }
    }
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
        [Vec3{ x: 0.0, y: -1.5, z: 0.0}, 
        Vec3{ x: 0.0, y: 1.5, z: 0.0}, 
        Vec3{ x: -2.0, y: 0.0, z: 0.0}, 
        Vec3{ x: 2.0, y: 0.0, z: 0.0}, 
        Vec3{ x: 0.0, y: 0.0, z: -2.0}, 
        Vec3{ x: 0.0, y: 0.0, z: 2.0}];

    let normals = [
        Vec3{ x: 0.0, y: 1.0, z: 0.0}, 
        Vec3{ x: 0.0, y: -1.0, z: 0.0}, 
        Vec3{ x: 1.0, y: 0.0, z: 0.0}, 
        Vec3{ x: -1.0, y: 0.0, z: 0.0}, 
        Vec3{ x: 0.0, y: 0.0, z: 1.0}, 
        Vec3{ x: 0.0, y: 0.0, z: -1.0}];

    for i in 0..points.len(){
        commands.spawn((components::PlaneBundle {
            point: components::Point(points[i]),
            normal: components::Normal(normals[i]),
        },
        components::Plane
        ));
    }

    commands.spawn(Camera3dBundle {
        transform: Transform::from_xyz(6.0, 2.0, -6.0).looking_at(Vec3::ZERO, Vec3::Y),
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
            radius: 0.1,
            subdivisions: 3,
        })
        .unwrap(),
    );

    for i in 0..10{
        for j in 0..1{
            for k in 0..10{

                let position = Vec3::new(
                    0.0 + (i as f32 * 0.1) + (j as f32 * 0.06),
                    0.0 + (j as f32 * 0.1),
                    0.0 + (k as f32 * 0.1) + (j as f32 * 0.06),
                );

                commands.spawn((components::ParticleBundle {
                    pbr: PbrBundle {
                        transform: Transform {
                            translation: position,
                            ..default()
                        },
                        mesh: mesh.clone(),
                        material: materials.add(Color::rgb(0.4, 0.98, 1.0).into()),
                        ..default()
                    },
                    mass: components::Mass(1.0),
                    density: components::Density(0.0),
                    pressure: components::Pressure(0.0),
                    acceleration: components::Acceleration(Vec3::ZERO),
                    velocity: components::Velocity(Vec3::ZERO),
                    last_pos: components::LastPos(position),
                },
                components::Particle,
                ));
            }
        }
    }
    

}

fn update_pressure(
    params: Res<FluidParameters>,
    mut query: Query<(
    &components::Mass,
    &Transform,
    &mut components::Density,
    &mut components::Pressure),
    With<components::Particle>>) {
    
    let mut iter = query.iter_combinations_mut();

    let mut delta: Vec3;
    let mut distance_sq: f32;

    // Set the density of each particle
    while let Some([(components::Mass(m1), transform1, mut density1, _pressure1), 
        (components::Mass(m2), transform2, mut density2, _pressure2)]) =
        iter.fetch_next()
    {
        delta = transform2.translation - transform1.translation;
        distance_sq = delta.length_squared();

        if distance_sq > params.smoothing_radius * params.smoothing_radius{
            continue;
        }

        let poly6 = get_poly6_smoothing(distance_sq, params.smoothing_radius);
        density1.0 += *m1 * poly6;
        density2.0 += *m2 * poly6;
    }

    // Set the pressure of the particle
    for (_mass, _transform, mut density, mut pressure) in &mut query {
        if density.0 < params.reference_density{
            density.0 = params.reference_density;
        }
        
        pressure.0 = params.pressure_constant * (density.0 - params.reference_density);
    }
}

fn update_acceleration(
    params: Res<FluidParameters>,
    mut query: Query<(
    &components::Mass,
    &Transform,
    &components::Density,
    &components::Pressure,
    &mut components::Acceleration,
    &components::Velocity),
    With<components::Particle>>){

    let mut iter = query.iter_combinations_mut();

    let mut rij: Vec3;
    let mut rji: Vec3;
    let mut distance_sq: f32;
    let mut spiky: f32;

    let mut vjvi: Vec3;
    let mut vivj: Vec3;
    let mut visc: f32;

    // Set the density of each particle
    while let Some([(components::Mass(m1), transform1, density1, pressure1, mut accel1, velocity1), 
        (components::Mass(m2), transform2, density2, pressure2, mut accel2, velocity2)]) =
        iter.fetch_next()
    {
        rij = transform1.translation - transform2.translation;
        distance_sq = rij.length_squared();

        if distance_sq > params.smoothing_radius * params.smoothing_radius{
            continue;
        }

        rij = rij.normalize();

        // vector from j to i is same as i to j * -1
        rji = rij * -1.0;

        spiky = get_spiky_smoothing(distance_sq, params.smoothing_radius);

        // Acceleration
        if density1.0.abs() > params.near_zero && density2.0.abs() > params.near_zero && !rji.is_nan(){

            // Removed mass since all particles have the same mass
            accel1.0 += -((m1 * (pressure1.0 + pressure2.0) / (2.0 * density1.0 * density2.0)))
                * spiky * rij;
            accel2.0 += -((m2 * (pressure1.0 + pressure2.0) / (2.0 * density1.0 * density2.0)))
            * spiky * rji;
        }

        // Viscocity Acceleration
        vjvi = velocity2.0 - velocity1.0;
        vivj = vjvi * -1.0;
        visc = get_viscosity_smoothing(distance_sq, params.smoothing_radius);

        if density2.0 > params.near_zero { 
            accel1.0 += params.viscosity_coef * (m2/density1.0) * (m1 / density2.0) * (vjvi) * visc;

        }
        if density1.0 > params.near_zero {
            accel2.0 += params.viscosity_coef * (m1/density2.0) * (m2 / density1.0) * (vivj) * visc;
        }
    }

    // Add gravity and clamp
    for (_mass, _transform, _density, _pressure, mut accel, _velocity) in &mut query {
        accel.0 += Vec3{ x: 0.0, y: -9.8, z: 0.0};

        if accel.0.length() > params.max_acceleration {
            accel.0 = accel.0.normalize() * params.max_acceleration;
        }
    }

}

fn plane_repulsion(
    mut particles: Query<(&components::Mass, 
        &mut components::Acceleration,
        &mut Transform),
        With<components::Particle>>,
    mut planes: Query<(&components::Point, &components::Normal), With<components::Plane>>,
){
    let mut distance_sq: f32 = 0.0;
    let mut rij: Vec3;  

    for (components::Mass(_m), mut accel, transform) in &mut particles {
        for (point, normal) in &mut planes{
            rij = point.0 * normal.0.abs() - (transform.translation * normal.0.abs());
            distance_sq = (rij).length_squared();

            let damp = accel.0.dot(normal.0.abs()).abs();

            accel.0 += normal.0 * wall_exp(distance_sq, 0.5, 1.1 * damp);
        }
    }
}

// A faux exponential so that the acceleration repulsion is stronger closer to the wall.
fn wall_exp(x: f32, smoothing_radius: f32, scale: f32) -> f32{
    let a: f32 = 500000.0;
    scale * ((a.powf((-(1.0 / smoothing_radius) * x) + 1.0) - 1.0) / (a - 1.0))
}

fn plane_collision(
    mut particles: Query<(&mut components::Velocity, &mut Transform), With<components::Particle>>,
    mut planes: Query<(&components::Point, &components::Normal), With<components::Plane>>,
){
    for (mut velocity, mut transform) in &mut particles {
        for (point, normal) in &mut planes{

            // If we are inside of the bounds or outside, but moving back into bounds continue
            if (transform.translation - point.0).dot(normal.0) > 0.0{
                continue;
            }

            velocity.0 = -((0.2 * normal.0).dot(velocity.0)) * (normal.0);
            transform.translation = ((point.0 - transform.translation) * normal.0.abs()) + transform.translation;
                
        }
        
    }
}

///
/// Input. distance_sq: The distance between two given particles
/// Ouput. The poly6 smoothing kernel
/// 
fn get_poly6_smoothing(distance_sq: f32, smoothing_radius: f32) -> f32{
    if distance_sq > smoothing_radius * smoothing_radius {
        return 0.0;
    }

    let mut poly6 = (smoothing_radius * smoothing_radius) - distance_sq;

    poly6 = f32::powf(poly6, 4.0);
    poly6 = (315.0 / (64.0 * PI * f32::powf(smoothing_radius, 12.0))) * poly6;

    poly6
}

///
/// Input. distance_sq: The distance between two given particles
/// Ouput. The poly6 smoothing kernel
/// 
fn get_spiky_smoothing(distance_sq: f32, smoothing_radius: f32) -> f32{
    let mut spiky = f32::powf(smoothing_radius - f32::sqrt(distance_sq), 2.0);
    spiky = (-45.0 / (PI * f32::powf(smoothing_radius, 6.0))) * spiky;

    spiky
}

fn get_viscosity_smoothing(distance_sq: f32, smoothing_radius: f32) -> f32{
    (45.0 / (PI * f32::powf(smoothing_radius, 6.0))) * (smoothing_radius - f32::sqrt(distance_sq))
}

fn integrate(
    params: Res<FluidParameters>,
    mut query: Query<(&mut components::Acceleration, &mut components::Velocity, &mut Transform)>) {

    for (mut acceleration, mut velocity, mut transform) in &mut query {
        velocity.0 = velocity.0 + (acceleration.0 * params.delta_time as f32);
        if velocity.0.length() > params.max_velocity {
            velocity.0 = velocity.0.normalize() * params.max_velocity;
        }

        transform.translation = transform.translation + (velocity.0 * params.delta_time as f32);
        acceleration.0 = Vec3::ZERO;
    }
}

fn clean_particle(mut query: Query<(&mut components::Acceleration, &mut components::Velocity, &mut components::Density, &mut components::Pressure)>){
    for (mut _acceleration, mut _velocity, mut density, mut pressure) in &mut query {
        density.0 = 0.0;
        pressure.0 = 0.0;
    }
}