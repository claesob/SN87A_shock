#FC = ifort
#FC = g77
FC = gfortran
FFLAGS =    -g -fbacktrace -Wall -ffixed-line-length-none   -fcheck=all  -mcmodel=medium
SRCS = main.f timed_hydro_shock.f trsub.f rad.f coll_rec.f radsub.f \
	radsub2.f pop_ch.f popsub_ch.f popca.f rates_1.f \
	rates_2.f big.f cross.f cfit.f cont_emiss2.f chianti_data.f chianti_data_2022_v4.f rechyfit.f verner_fit.f coll_ionization_verner.f auger.f read_auger.f num_shell.f read_rr_diel_badnell_v2.f dielbadnell.f 
OBJS = $(SRCS:.f=.o)

shock : $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS)

clean:
	-rm $(OBJS)






