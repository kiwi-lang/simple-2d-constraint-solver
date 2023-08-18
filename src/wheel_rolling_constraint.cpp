#include "../include/wheel_rolling_constraint.h"

#include <cmath>

namespace atg_scs {
void WheelRollingConstraint::calculate(Output* output, atg_scs::SystemState *state) {

    const int rollingWheel = m_bodies[0]->index;
    const int slidingWheel = m_bodies[1]->index;

    const double q6 = state->theta[slidingWheel];

    const double q1_dot = state->v_x[rollingWheel];
    const double q2_dot = state->v_y[rollingWheel];
    const double q3_dot = state->v_theta[rollingWheel];

    const double q6_dot = state->v_theta[slidingWheel];

    const double cos_q6 = std::cos(q6);
    const double sin_q6 = std::sin(q6);

    const double u_x = cos_q6;
    const double u_y = sin_q6;

    const double u_x_dot = -sin_q6 * q6_dot;
    const double u_y_dot = cos_q6 * q6_dot;

    const double contactVelocity_x = q3_dot * m_wheelRadius * u_x + q1_dot;
    const double contactVelocity_y = q3_dot * m_wheelRadius * u_y + q2_dot;
    const double contactVelocity_2 = contactVelocity_x * contactVelocity_x  + contactVelocity_y * contactVelocity_y;

    const double friction = contactVelocity_2 > 0.5 ? m_dynamicFriction: m_staticFriction;

    output->J[0][0] = u_x;
    output->J[0][1] = u_y;
    output->J[0][2] = m_wheelRadius;

    output->J_dot[0][0] = u_x_dot;
    output->J_dot[0][1] = u_y_dot;
    output->J_dot[0][2] = 0;

    output->C[0] = 0;

    output->v_bias[0] = 0;

    output->limits[0][0] = -friction * m_wheelRadius;
    output->limits[0][1] = friction * m_wheelRadius;

    output->kd[0] = m_kd;
    output->ks[0] = m_ks;
}
}