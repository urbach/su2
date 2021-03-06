strcpl <- function(beta, d) {
  2.*(beta/4. 
    - (beta^3)/32. + 6.*(d/192.-11./1152.)*beta^5
    + 8.*(-d/256. + 757./98304.)*beta^7
    + 10.*(d^2/512. - d*85./12288. + 2473./409600.)*beta^9
    + 12.*(-d^2*29./12288. + d*2467./262144. - 1992533./212336640.)*beta^11
    + 14.*((d^3)*5./4096. - d^2*237./32768. + d*178003./11796480.-38197099./3468165120.)*beta^13
    + 16.*(-(d^3)*15./8192. + d^2*1485./131072. -d*53956913./2264924160. + 11483169709./676457349120.)*beta^15 )
}

weakcpl <- function(beta, d, L=2) {
  E = 1.-1./L**d

  d0=E/d
  d1=0.5*(E^2)/(d^2)
  d2=(2*E^3)/(3.*d^3) + 0.2442 /(3.*2^4)
  d3=(9.*E^4)/(8.*d^4) + 0.1859/(2.^6) + 0.2442*E/(16.*d*d*(d-1.)) 
  
  return(1.-d0/beta -d1/beta^2 - d2/beta^3 - d3/beta^4)
}
