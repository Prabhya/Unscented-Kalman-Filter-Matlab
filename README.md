# Unscented-Kalman-Filter-Matlab
Unscented kalman filter uses the second order approximation of the taylor series expansion for nonlinear systems.

The unscented filter can work properly even if frequency of the observation data provided is not high enough as compared to EKF.
The 2 reason for this being UKF uses the second order approximation compared to ekf which uses only first order approximation of taylor series.
Second reason being augmentation of states, this being, the noise and unknown states can all be augmented into one state called the augmented state. This helps in propagating the state keeping noise into higher consideration.
