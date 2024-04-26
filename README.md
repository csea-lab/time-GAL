![Time-GAL toolbox Banner](pictures/TGAL_banner.png)
![Visits Badge](https://badges.pufler.dev/visits/csea-lab/time-GAL?style=for-the-badge)
<img alt="Version" src="https://img.shields.io/badge/Version-1.0.0-blue?style=for-the-badge">
<img alt="Version" src="https://img.shields.io/badge/Language-MATLAB-orange?style=for-the-badge">
<br>
<br>
The Time-GAL toolbox offers a framework to investigate how a  cognitive  processing is shared within different brain regions. This toolbox computes the Generalization Across Location (GAL), a decoding-based methodology applying multi-variate pattern analysis (MVPA) to discriminate decodable information inlaid within the neural activity of the brain. Instead of using spatial information from electrodes or voxels, this methodology leverages the temporal neural dynamics, i.e. recording time-series, using this dimension as features for the MVPA-ased classifier model.

Conceptually, if a time-based decoder trained at location (channel, sensor or voxel) also classifies the data at other location, then we interpret this as evidence that the temporal dynamics at both locations contain shared decodable information about the conditions—a non-parametric way of describing spatial interdependency or connectivity. Therefore, we can visualize how different brain areas interacts between them during a cognitive or affective related task.

Furthermore, this toolbox computes not only the avobe mentioned spatial decoding backward model, but also a feedforward model of the temporal neural activity during the trial. Combining both, we can inspect how neural dynamics elaps over the trial showing not only the topographical brain aeras involved in the cognitive process, but also the precise generalization connectivity patterns between the brain areas across time.

# How it works

In the decades, 
