function Y = correct_instrument_response(XL, Y)
% XL - x vector in nm
% Y - Intensity values to correct

load DilorXY_Instrument_Response.mat;

Y = Y./interp1(instrument_response(:,1),instrument_response(:,2),XL);


