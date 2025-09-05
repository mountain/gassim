#![allow(missing_docs)] // Doc comments are added on public items in submodules per AGENT.md

//! Core simulation data structures and foundational types for GasSim.
//!
//! This module defines the foundational types needed by the event-driven
//! molecular dynamics engine, aligning with Phase 0 of docs/plan.md.

pub mod event;
pub mod particle;

pub use event::{Event, EventKind};
pub use particle::Particle;
