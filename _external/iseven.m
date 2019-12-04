%> @file iseven.m
%> @brief  Returns true if the input is an even number. Syntax: tf = iseven(x)
%> @author Dan Kominsky, Copyright 2012  Prime Photonics, LC.
%> @date   2012
%======================================================================
%> @param tf = The number of modes
%> @retval x	 = numeric (real) input of any number of dimensions and sizes
%> @retval Mj	 =  true for each element that is an even whole number
% ======================================================================
  function tf = iseven(x)
  if ~isreal(x)
    error('iseven:badinput','iseven requires real inputs');
  else
    tf = mod(x,2)==0;
  end
end

  % iseven - Returns true if the input is an even whole number
  % Syntax: tf = iseven(x)
  % x - numeric (real) input of any number of dimensions and sizes
  % tf - true for each element that is an even whole number
  %   Example
  %  tf = iseven([1,3,4]);
  %
  %   See also: isodd
  
  % AUTHOR    : Dan Kominsky
  % Copyright 2012  Prime Photonics, LC.
