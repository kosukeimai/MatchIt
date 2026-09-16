// // [[Rcpp::depends(RcppProgress)]]
// #include <RcppArmadillo.h>
// #include <progress.hpp>
// // [[Rcpp::depends(RcppArmadillo)]]
// // using namespace Rcpp;
//
// // [[Rcpp::export]]
// arma::uvec energy_match_attC(const arma::vec& t,
//                        const arma::mat d,
//                        const bool& disl_prog = false) {
//
//   //Energy distance formula:
//   //2/n1n0 sum(D[t,c]) - 1/n1n1 sum(D[t,t]) - 1/n0n0 sum(D[c,c])
//
//   arma::uvec treated_ind = arma::find(t == 1);
//   arma::uvec control_ind = arma::find(t == 0);
//
//   double n1 = treated_ind.size();
//   double n0 = control_ind.size(), n0_ = control_ind.size();
//
//   // progress bar
//   int prog_length = n0;
//   Progress p(prog_length, disl_prog);
//
//   //Compute each unit's contribution to energy dist
//   //Find unit that contributes maximally to energy distance
//   //Remove that unit
//
//   arma::uvec drop_order(n0 + 1);
//   drop_order.fill(NA_INTEGER);
//   arma::vec Ys(n0 + 1);
//   Ys.fill(NA_REAL);
//
//   arma::mat d10 = d(treated_ind, control_ind);
//   arma::vec cSd10 = arma::sum(d10, 0);
//   double Sd10 = arma::accu(cSd10);
//
//   arma::mat d11 = d(treated_ind, treated_ind);
//   double Sd11 = arma::accu(d11);
//
//   arma::mat d00 = d(control_ind, control_ind);
//   arma::vec cSd00 = arma::sum(d00, 0);
//   double Sd00 = arma::accu(cSd00);
//   arma::uvec d00i = arma::conv_to<arma::uvec>::from(arma::linspace(0, n0 - 1, n0));
//
//   Ys[0] = 2*Sd10/(n1*n0) - Sd11/(n1*n1) - Sd00/(n0*n0);
//
//   double min_edist = Ys[0];
//   double edist;
//
//   arma::vec control_contributions =
//     (2/(n0*n1) - 2/((n0-1)*n1))*Sd10 +
//     (-1/pow(n0,2) + 1/pow((n0-1),2))*Sd00 +
//     (2/((n0-1)*n1))*cSd10 -
//     (2/pow((n0-1),2))*cSd00;
//
//   // double largest_contribution(n0);
//   int drop_control, control_ind_to_drop;
//   arma::uvec d00i_control_dropped(n0), d00i_control_kept(n0), i0(n0);
//
//   // double cc = control_contributions;
//   // Rcout << "control_contributions = " << cc << "\n";
//
//   int k = 1;
//   while(true) {
//
//     i0 = arma::conv_to<arma::uvec>::from(arma::linspace(0, n0 - 1, n0));
//
//     // largest_contribution = max(control_contributions);
//     //Find which units have the largest contribution to the energy distance
//     drop_control = control_contributions.index_max();
//
//     control_ind_to_drop = control_ind[drop_control];
//
//     //Remove those units and decrease
//     control_ind = control_ind.elem(arma::find(i0 != drop_control));
//     // control_ind.erase(drop_control);
//     n0 = control_ind.size();
//
//     //If removed units are last units, don't remove and stop
//     if (n0 - 1 <= 0) break;
//
//     p.update(n0_ - n0);
//
//     //Compute new edist with dropped units removed
//
//     d00i_control_dropped = d00i[drop_control];
//     d00i = d00i.elem(arma::find(i0 != drop_control));
//     // d00i.erase(drop_control);
//
//     //d10 contributions
//     Sd10 = Sd10 - arma::accu(d10.cols(d00i_control_dropped));
//     cSd10 = cSd10.elem(arma::find(i0 != drop_control));
//     // cSd10.erase(drop_control);
//
//     //d00 contributions
//     Sd00 = Sd00 - 2*arma::accu(d00(d00i, d00i_control_dropped)) -
//       arma::accu(d00(d00i_control_dropped, d00i_control_dropped));
//
//     cSd00 = cSd00.elem(arma::find(i0 != drop_control));
//     // cSd00.erase(drop_control);
//     cSd00 = cSd00 - arma::sum(d00(d00i_control_dropped, d00i), 0);
//     // d00i = d00i_control_kept;
//
//     edist = 2*Sd10/(n1*n0) - Sd11/(n1*n1) - Sd00/(n0*n0);
//
//     //After passing sample size threshold of 90% of original N, stop if
//     //new edist is larger than smallest edist
//     if (edist < min_edist) min_edist = edist;
//     // else if (n0/n0_ < .9 && edist - min_edist > .2*(Ys[0] - min_edist)) break;
//
//     //Record new edist and units dropped
//     Ys[k] = edist;
//     drop_order[k] = control_ind_to_drop;
//
//     //Compute unit contributions after having dropped units
//     control_contributions =
//       (2/(n0*n1) - 2/((n0-1)*n1))*Sd10 +
//       (-1/pow(n0,2) + 1/pow(n0-1,2))*Sd00 +
//       (2/((n0-1)*n1))*cSd10 -
//       (2/pow(n0-1, 2))*cSd00;
//
//     k++;
//   }
//
//   // Ys = Ys[seq(0, k - 2)];
//   // drop_order = drop_order[seq(0, k - 2)];
//
//   p.update(prog_length);
//
//   return drop_order + 1;
// }