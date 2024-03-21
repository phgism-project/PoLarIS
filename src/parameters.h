#ifndef PARAMETERS_
/* Constants */

/*	 Paterson & budd (1982) */
#  define	 ARRHENIUS_T	263.15
#  define	 ARRHENIUS_A0   3.61E-13     
#  define	 ARRHENIUS_A1	1.73E3
#  define	 ARRHENIUS_Q0	6.0E4
#  define	 ARRHENIUS_Q1	13.9E4

#  define	 RHO_ICE		910	/* Density of ice */
#  define	 RHO_WATER		1000	/* Density of sea water */
#  define	 GRAVITY		9.81	/* Acceleration due to gravity */
#  define	 POWER_N		3	/* power in Glen's law */
#  define	 TEMP_WATER		273.15	/* Triple point of water */
#  define	 GEOTHE_FLUX		3E-2	/* Geothermal heat flux */
#  define	 THEM_CONDUCT		2.1	/* Thermal conductivity of ice */
#  define	 HEAT_CAPACITY		2009	/* Heat capacity of ice  */
#  define	 BETA_MELT		8.66E-4	/* Dependence of melting on depth */
#  define	 LATENT_CAPACITY	3.35E5	/* Latent-heat capacity of ice */
#  define	 GAS_CONSTANT		8.314	/* Gas constat */
#  define	 SEC_PER_YEAR		31556926	/* Seconds per year */
#  define        MIN_EFFECTIVE_STRAIN   1e-15	/* Effective strian minimum */
#  define        HEAT_SOURCE            0       /* Volume heat soruce */

#  define 	 HEIGHT_EPS     (1./LEN_SCALING)	/* Height eps 1 meter */


#define PARAMETERS_
#endif
