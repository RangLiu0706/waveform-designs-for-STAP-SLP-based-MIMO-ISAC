### About the paper
This is a code package for the paper: 
R. Liu, M. Li, Q. Liu, and A. L. Swindlehurst, "Joint waveform and filter designs for STAP-SLP-based MIMO-DFRC systems,”IEEE J. Sel. Areas Commun., vol. 40, no. 6, pp. 1918-1931, Jun. 2022.

@ARTICLE{9724259,
  author={Liu, Rang and Li, Ming and Liu, Qian and Swindlehurst, A. Lee},
  journal={IEEE Journal on Selected Areas in Communications}, 
  title={Joint Waveform and Filter Designs for STAP-SLP-Based MIMO-DFRC Systems}, 
  year={2022},
  volume={40},
  number={6},
  pages={1918-1931},
  keywords={Radar;Clutter;Sensors;Radar clutter;Quality of service;Signal to noise ratio;Wireless communication;Dual-functional radar-communication (DFRC);integrated sensing and communication (ISAC);space-time adaptive processing (STAP);symbol-level precoding (SLP);multi-input multi-output (MIMO)},
  doi={10.1109/JSAC.2022.3155501}}


- If you use this simulation code package in any way, please cite the original paper above.
- All codes are contributed by Rang Liu (email: rangl2@uci.edu; website: https://rangliu0706.github.io/). 
   Please feel free to contact with her if you have any suggestions. 
- The link of this paper is: https://ieeexplore.ieee.org/document/9769997
- More information can be found at: https://www.minglabdut.com/resource.html
- Copyright Notice: This code is licensed for personal, non-commercial use only, specifically for academic purposes. Copyright reserved by the MingLab (led by Prof. Ming Li), School of Information and Communication Engineering, Dalian University of Technology, Dalian 116024, China. 


### Software platform
- Please note that the MATLAB2022b is used for this simulation code package, and there may be some imcompatibility problems among different sofrware versions. 
- To run those codes, please download and install [CVX](http://cvxr.com/cvx/) & [Manopt](https://www.manopt.org/)

### Content of this simulation code package
- The files "main_iterations", "main_SINR_SNR", "main_SINR_Ku", "main_SINR_P", "main_SINR_Nt", "main_SINR_fd", "main_SINR_ratio", and "main_SINR_epsi" are used to generate Figs. 3-10, respectively.

Abstract of the paper: 
Dual-function radar-communication (DFRC), which can simultaneously perform both radar and communication functionalities using the same hardware platform, spectral resource and transmit waveform, is a promising technique for realizing integrated sensing and communication (ISAC). Space-time adaptive processing (STAP) in multi-antenna radar systems is the primary tool for detecting moving targets in the presence of strong clutter. The idea of joint spatial-temporal optimization in STAP-based radar systems is consistent with the concept of symbol-level precoding (SLP) for multi-input multi-output (MIMO) communications, which optimizes the transmit waveform for each of the transmitted symbols. In this paper, we combine STAP and SLP and propose a novel STAP-SLP-based DFRC system that enjoys the advantages of both techniques. The radar output signal-to-interference-plus-noise ratio (SINR) is maximized by jointly optimizing the transmit waveform and receive filter, while satisfying the communication quality-of-service (QoS) constraint and various waveform constraints including constant-modulus, similarity and peak-to-average power ratio (PAPR). An efficient algorithm framework based on majorization-minimization (MM) and nonlinear equality constrained alternative direction method of multipliers (neADMM) methods is proposed to solve these complicated non-convex optimization problems. Simulation results verify the effectiveness of the proposed STAP-SLP-based MIMO-DRFC scheme and the associate algorithms.

