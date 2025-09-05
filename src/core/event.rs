use crate::error::{Error, Result};
use ordered_float::NotNan;
use std::cmp::Ordering;

/// Kinds of events that can occur in the EDMD engine.
///
/// Tie-breaking for deterministic ordering prefers `P2P` < `P2W` when times are equal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventKind {
    /// Particle-to-particle collision between particles `i` and `j`.
    P2P { i: u32, j: u32 },
    /// Particle-to-wall collision for particle `i` and wall `wall_id`.
    P2W { i: u32, wall_id: u32 },
}

impl EventKind {
    #[inline]
    fn order_key(&self) -> (u8, u32, u32) {
        match *self {
            EventKind::P2P { i, j } => (0, i, j),
            EventKind::P2W { i, wall_id } => (1, i, wall_id),
        }
    }
}

/// A scheduled event in the priority queue with deterministic ordering.
///
/// - `time`: event occurrence time (finite, non-NaN).
/// - `kind`: event kind and participants.
/// - `cc_i`, `cc_j`: collision-count snapshots for invalidation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Event {
    pub time: NotNan<f64>,
    pub kind: EventKind,
    pub cc_i: u64,
    pub cc_j: Option<u64>,
}

impl Event {
    /// Create a new event, validating that time is finite and non-NaN.
    pub fn new(time: f64, kind: EventKind, cc_i: u64, cc_j: Option<u64>) -> Result<Self> {
        if time.is_nan() {
            return Err(Error::InvalidParam("event time cannot be NaN".into()));
        }
        if !time.is_finite() {
            return Err(Error::InvalidParam("event time must be finite".into()));
        }
        let time = NotNan::new(time)
            .map_err(|_| Error::InvalidParam("event time cannot be NaN".into()))?;
        Ok(Self {
            time,
            kind,
            cc_i,
            cc_j,
        })
    }

    /// Returns the raw f64 event time.
    #[inline]
    pub fn time_f64(&self) -> f64 {
        self.time.into_inner()
    }

    /// Validate against current collision-counts. If the second participant does not exist
    /// (e.g., a P2W event), pass `None` for `cc_j_now`.
    #[inline]
    pub fn is_valid(&self, cc_i_now: u64, cc_j_now: Option<u64>) -> bool {
        if self.cc_i != cc_i_now {
            return false;
        }
        match (self.cc_j, cc_j_now) {
            (Some(a), Some(b)) => a == b,
            (None, None) => true,
            // Event with only `i` remains valid regardless of unrelated `cc_j_now`.
            (None, Some(_)) => true,
            // If the event expected a second participant but none is provided now, invalidate.
            (Some(_), None) => false,
        }
    }
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.time.cmp(&other.time) {
            Ordering::Equal => {
                let a = self.kind.order_key();
                let b = other.kind.order_key();
                match a.cmp(&b) {
                    Ordering::Equal => {
                        // Final tie-breaker on cc snapshots to ensure a total order.
                        (self.cc_i, self.cc_j.unwrap_or(0))
                            .cmp(&(other.cc_i, other.cc_j.unwrap_or(0)))
                    }
                    o => o,
                }
            }
            o => o,
        }
    }
}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use EventKind::{P2P, P2W};

    #[test]
    fn new_event_rejects_nan_time() {
        let err = Event::new(f64::NAN, P2P { i: 1, j: 2 }, 0, Some(0)).unwrap_err();
        assert!(err.to_string().contains("NaN"));
    }

    #[test]
    fn ordering_by_time() -> Result<()> {
        let e1 = Event::new(1.0, P2P { i: 0, j: 1 }, 0, Some(0))?;
        let e2 = Event::new(2.0, P2W { i: 0, wall_id: 0 }, 0, None)?;
        assert!(e1 < e2);
        Ok(())
    }

    #[test]
    fn tie_breaker_prefers_p2p_over_p2w() -> Result<()> {
        let t = 5.0;
        let a = Event::new(t, P2P { i: 0, j: 1 }, 3, Some(4))?;
        let b = Event::new(t, P2W { i: 0, wall_id: 1 }, 3, None)?;
        assert!(a < b); // P2P < P2W at equal time
        Ok(())
    }

    #[test]
    fn is_valid_checks_collision_counts() -> Result<()> {
        let e_pp = Event::new(1.0, P2P { i: 1, j: 2 }, 10, Some(20))?;
        assert!(e_pp.is_valid(10, Some(20)));
        assert!(!e_pp.is_valid(11, Some(20)));
        assert!(!e_pp.is_valid(10, Some(21)));
        assert!(!e_pp.is_valid(10, None));

        let e_pw = Event::new(1.0, P2W { i: 3, wall_id: 0 }, 7, None)?;
        assert!(e_pw.is_valid(7, None));
        // Passing an unrelated second cc should not invalidate a P2W event.
        assert!(e_pw.is_valid(7, Some(999)));
        Ok(())
    }
}
