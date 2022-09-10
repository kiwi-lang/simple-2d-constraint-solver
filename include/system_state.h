#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SYSTEM_STATE_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SYSTEM_STATE_H

#include "wrapped_pointer.h"

namespace atg_scs {
    class SystemState {
        public:
            SystemState();
            ~SystemState();

            void copy(const SystemState *state);
            void resize(int bodyCount, int constraintCount);
            void destroy();

            void localToWorld(double x, double y, double *x_t, double *y_t, int body);
            void velocityAtPoint(double x, double y, double *v_x, double *v_y, int body);
            void applyForce(double x_l, double y_l, double f_x, double f_y, int body);

            Ptr<int> indexMap;

            Ptr<double> a_theta;
            Ptr<double> v_theta;
            Ptr<double> theta;

            Ptr<double> a_x;
            Ptr<double> a_y;
            Ptr<double> v_x;
            Ptr<double> v_y;
            Ptr<double> p_x;
            Ptr<double> p_y;

            Ptr<double> f_x;
            Ptr<double> f_y;
            Ptr<double> t;

            Ptr<double> r_x;
            Ptr<double> r_y;
            Ptr<double> r_t;

            Ptr<double> m;

            int n;
            int n_c;
            double dt;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SYSTEM_STATE_H */
