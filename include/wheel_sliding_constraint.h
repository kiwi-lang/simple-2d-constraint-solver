
#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_WHEELSLIDINGCONSTRAINT_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_WHEELSLIDINGCONSTRAINT_H

#include "constraint.h"

namespace atg_scs {
    class WheelSlidingConstraint : public Constraint {
        public:
            WheelSlidingConstraint();
            virtual ~WheelSlidingConstraint();
            
            virtual void calculate(Output *output, SystemState *system);

            double m_kd;
            double m_ks;
            double m_staticFriction;
            double m_dynamicFriction;
            double m_wheelRadius;
            RigidBody* m_rollingWheel;

    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_LINK_CONSTRAINT_H */
