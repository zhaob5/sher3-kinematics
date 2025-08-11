# SHER-3.0 Kinematics and Jacobian Analysis

This repository contains detailed derivations and explanations of the **forward kinematics**, **inverse kinematics**, **Jacobian**, and **Jacobian Inverse** analysis for the **SHER-3.0 surgical robot**. It includes analytical computation results in Python, validation through SolidWorks motion simulation, and detailed documentation for future insights and reference.

![](https://github.com/zhaob5/sher3-kinematics/blob/main/Figures/SHER-3.0.png)

## Recommended reading order:

1. **[Kinematics of a 3-PUU Delta-Style Parallel Robot](https://github.com/zhaob5/sher3-kinematics/blob/main/Kinematics%20of%20a%203-PUU%20Delta-Style%20Parallel%20Robot.pdf)** – foundational understanding of the delta platform.

2. **[Kinematics of the Roll-Tilt Mechanism in SHER-3.0](https://github.com/zhaob5/sher3-kinematics/blob/main/Kinematics%20of%20the%20Roll-Tilt%20Mechanism%20in%20SHER-3.0.pdf)** – breakdown of the rotating four-bar linkage structure.
   
3. **[(Optional) SOLIDWORKS Motion Study Tutorial: Setup and Kinematic Analysis for SHER-3.0](https://github.com/zhaob5/sher3-kinematics/blob/main/Kinematics%20of%20the%20Roll-Tilt%20Mechanism%20in%20SHER-3.0.pdf)** - a quick-start reference for running motion simulation in SOLIDWORKS.

4. **[Kinematics and Jacobian Analysis of SHER-3.0](https://github.com/zhaob5/sher3-kinematics/blob/main/Kinematics%20and%20Jacobian%20Analysis%20of%20SHER-3.0.pdf)** – full system modeling and differential kinematics.

5. **[Closed-Form Jacobian Inverse Derivation for SHER-3.0](https://github.com/zhaob5/sher3-kinematics/blob/main/Closed-Form%20Jacobian%20Inverse%20Derivation%20for%20SHER-3.0.pdf)** - derivation of the Jacobian Inverse with physical intuition.

## License

Copyright © 2025 JHU LCSR. Developed by [AMIRo Research Lab](https://amiro.lcsr.jhu.edu/).

See the [LICENSE](./LICENSE) file for details.

## Disclaimer

The derivation and analysis presented in this repository are based on analytical methods applied to a custom-designed robotic system. Readers are strongly encouraged to focus on the **methodology and reasoning** behind the derivations rather than relying solely on the final equations.

**Please independently verify all results before applying them to any critical or safety-sensitive applications.** This work is provided “as is” without warranty of any kind.
