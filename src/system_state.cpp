#include "../include/system_state.h"

#include "../include/utilities.h"

#include <assert.h>
#include <cstring>
#include <cmath>

atg_scs::SystemState::SystemState() {
    indexMap = nullptr;

    a_theta = nullptr;
    v_theta = nullptr;
    theta = nullptr;

    a_x = nullptr;
    a_y = nullptr;
    v_x = nullptr;
    v_y = nullptr;
    p_x = nullptr;
    p_y = nullptr;

    f_x = nullptr;
    f_y = nullptr;
    t = nullptr;

    m = nullptr;

    r_x = 0;
    r_y = 0;
    r_t = 0;

    n = 0;
    n_c = 0;
    dt = 0.0;
}

atg_scs::SystemState::~SystemState() {
    assert(n == 0);
    assert(n_c == 0);
}

void atg_scs::SystemState::copy(const SystemState *state) {
    resize(state->n, state->n_c);

    if (state->n == 0) {
        return;
    }

    std::memcpy((void *)indexMap.get(), (void*)state->indexMap.get(), sizeof(int) * n_c);

    std::memcpy((void *)a_theta.get(), (void *)state->a_theta.get(), sizeof(double) * n);
    std::memcpy((void *)v_theta.get(), (void *)state->v_theta.get(), sizeof(double) * n);
    std::memcpy((void *)theta.get(), (void *)state->theta.get(), sizeof(double) * n);

    std::memcpy((void *)a_x.get(), (void *)state->a_x.get(), sizeof(double) * n);
    std::memcpy((void *)a_y.get(), (void *)state->a_y.get(), sizeof(double) * n);
    std::memcpy((void *)v_x.get(), (void *)state->v_x.get(), sizeof(double) * n);
    std::memcpy((void *)v_y.get(), (void *)state->v_y.get(), sizeof(double) * n);
    std::memcpy((void *)p_x.get(), (void *)state->p_x.get(), sizeof(double) * n);
    std::memcpy((void *)p_y.get(), (void *)state->p_y.get(), sizeof(double) * n);

    std::memcpy((void *)f_x.get(), (void *)state->f_x.get(), sizeof(double) * n);
    std::memcpy((void *)f_y.get(), (void *)state->f_y.get(), sizeof(double) * n);
    std::memcpy((void *)t.get(), (void *)state->t.get(), sizeof(double) * n);

    std::memcpy((void *)m.get(), (void *)state->m.get(), sizeof(double) * n);

    std::memcpy((void *)r_x.get(), (void *)state->r_x.get(), sizeof(double) * n_c * 2);
    std::memcpy((void *)r_y.get(), (void *)state->r_y.get(), sizeof(double) * n_c * 2);
    std::memcpy((void *)r_t.get(), (void *)state->r_t.get(), sizeof(double) * n_c * 2);
}

void atg_scs::SystemState::resize(int bodyCount, int constraintCount) {
    assert(bodyCount > 0 && constraintCount > 0);

    if (n >= bodyCount && n_c >= constraintCount) {
        return;
    }

    destroy();

    n = bodyCount;
    n_c = constraintCount;

    indexMap.make(n_c);

    a_theta.make(n);
    v_theta.make(n);
    theta.make(n);

    a_x.make(n);
    a_y.make(n);
    v_x.make(n);
    v_y.make(n);
    p_x.make(n);
    p_y.make(n);

    f_x.make(n);
    f_y.make(n);
    t.make(n);

    m.make(n);

    r_x.make(n_c * 2);
    r_y.make(n_c * 2);
    r_t.make(n_c * 2);
}

void atg_scs::SystemState::destroy() {
    if (n > 0) {
        (a_theta.destroy());
        (v_theta.destroy());
        (theta.destroy());

        (a_x.destroy());
        (a_y.destroy());
        (v_x.destroy());
        (v_y.destroy());
        (p_x.destroy());
        (p_y.destroy());

        (f_x.destroy());
        (f_y.destroy());
        (t.destroy());

        (m.destroy());
    }

    if (n_c > 0) {
        (indexMap.destroy());

        (r_x.destroy());
        (r_y.destroy());
        (r_t.destroy());
    }

    n = 0;
    n_c = 0;
}

void atg_scs::SystemState::localToWorld(
        double x,
        double y,
        double *x_t,
        double *y_t,
        int body)
{
    const double x0 = p_x[body];
    const double y0 = p_y[body];
    const double theta = this->theta[body];

    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);

    *x_t = cos_theta * x - sin_theta * y + x0;
    *y_t = sin_theta * x + cos_theta * y + y0;
}

void atg_scs::SystemState::velocityAtPoint(
        double x,
        double y,
        double *v_x,
        double *v_y,
        int body)
{
    double w_x, w_y;
    localToWorld(x, y, &w_x, &w_y, body);

    const double v_theta = this->v_theta[body];
    const double angularToLinear_x = -v_theta * (w_y - this->p_y[body]);
    const double angularToLinear_y = v_theta * (w_x - this->p_x[body]);

    *v_x = this->v_x[body] + angularToLinear_x;
    *v_y = this->v_y[body] + angularToLinear_y;
}

void atg_scs::SystemState::applyForce(
    double x_l,
    double y_l,
    double f_x,
    double f_y,
    int body)
{
    double w_x, w_y;
    localToWorld(x_l, y_l, &w_x, &w_y, body);

    this->f_x[body] += f_x;
    this->f_y[body] += f_y;

    this->t[body] +=
        (w_y - this->p_y[body]) * -f_x +
        (w_x - this->p_x[body]) * f_y;
}
