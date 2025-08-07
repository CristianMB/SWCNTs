% dipbackground  Intermolecular dipolar background decay
%
%  Vinter = dipbackground(t,conc,lambda)
%  dipbackground(...)
%
%  Calculates the intermolecular dipolar background decay over the time axis,
%  for spin concentration conc (in µM) and modulation depth lambda (between 0 and 1).
%
%  Input:
%    t       time axis (µs)
%    conc    spin concentration (µM)
%    lambda  modulation depth
%
%  Output:
%    Vinter  intermolecular background decay function
%
%  Example:
%    t = linspace(-0.2,5);  % µs
%    Vinter = dipbackground(t,100,0.3);
%    plot(t,Vinter);
