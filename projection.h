#pragma once

int my_lrs_init();
// void lrs_print_data(lrs_dic *P);
struct GMPmat *projection(struct GMPmat *inp, int d);
struct GMPmat *H2V(struct GMPmat *inp);
struct GMPmat *V2H(struct GMPmat *inp);
struct GMPmat *reducemat(struct GMPmat *inp);
struct GMPmat *reducevertices(struct GMPmat *inp);