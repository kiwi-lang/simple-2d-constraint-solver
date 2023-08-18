#include "../include/wheel_sliding_constraint.h"

#include <cmath>

namespace atg_scs {
void WheelSlidingConstraint::calculate(Output* output, atg_scs::SystemState *state) {

    const int slidingWheel = m_bodies[0]->index;
    const int rollingWheel = m_rollingWheel->index;

    const double q1 = state->p_x[slidingWheel];
    const double q2 = state->p_y[slidingWheel];
    const double q3 = state->theta[slidingWheel];

    const double q1_dot = state->v_x[slidingWheel];
    const double q2_dot = state->v_y[slidingWheel];
    const double q3_dot = state->v_theta[slidingWheel];

    const double rollingVelocity = state->v_theta[rollingWheel];

    const double cos_q3 = std::cos(q3);
    const double sin_q3 = std::sin(q3);

    const double u_x = cos_q3;
    const double u_y = sin_q3;

    const double u_prep_x = -sin_q3;
    const double u_prep_y = cos_q3;

    const double u_prep_x_dot = -cos_q3 * q3_dot;
    const double u_prep_y_dot = -sin_q3 * q3_dot;

    const double contactVelocity_x = rollingVelocity * m_wheelRadius * u_x + q1_dot;
    const double contactVelocity_y = rollingVelocity * m_wheelRadius * u_y + q2_dot;
    const double contactVelocity_2 = contactVelocity_x * contactVelocity_x + contactVelocity_y * contactVelocity_y;

    const double friction = contactVelocity_2 > 0.5 ? m_dynamicFriction: m_staticFriction;

    output->J[0][0] = u_prep_x;
    output->J[0][1] = u_prep_y;
    output->J[0][2] = 0;

    output->J_dot[0][0] = u_prep_x_dot;
    output->J_dot[0][1] = u_prep_y_dot;
    output->J_dot[0][2] = 0;

    output->C[0] = 0;

    output->v_bias[0] = 0;

    output->limits[0][0] = -friction;
    output->limits[0][1] = friction;

    output->kd[0] = m_kd;
    output->ks[0] = m_ks;
}
}