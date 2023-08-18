
#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_WheelRollingConstraint_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_WheelRollingConstraint_H

#include "constraint.h"

namespace atg_scs {
    class WheelRollingConstraint : public Constraint {
        public:
            WheelRollingConstraint();
            virtual ~WheelRollingConstraint();
            
            virtual void calculate(Output *output, SystemState *system);

            double m_kd;
            double m_ks;
            double m_staticFriction;
            double m_dynamicFriction;
            double m_wheelRadius;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_LINK_CONSTRAINT_H */
