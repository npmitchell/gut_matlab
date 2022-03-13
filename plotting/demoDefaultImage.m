function demoDefaultImage()

    disp('loading default image')
    defImage = pow2(get(0,'DefaultImageCData'),47);
    r = bitslice(defImage,0,0);
    g = bitslice(defImage,17,17);
    b = bitslice(defImage,34,34);
    imshow(cat(3,r,g,b));
    pause(1)
    
    % Make montage
    disp('prepping montage')
    imgCell = repmat({zeros(size(defImage))},8,7);
    for shift = 0:52
      imgCell{shift+1} = bitrotate(defImage,shift);
      imagesc(bitrotate(defImage,shift)) ;
      colormap gray
      pause(0.5)
    end
    allImages = cell2mat(imgCell.');
    imshow(allImages,[min(allImages(:)) max(allImages(:))]);
end

function b = bitslice(a,lowbit,highbit)
% BITSLICE(A,LOWBIT,HIGHBIT)
%
%   Scales return values so that the maximum is 1.0.

b = a / 2^(lowbit);
b = fix(b);
b = b / 2^(highbit - lowbit + 1);
b = b - fix(b);

b = b / max(b(:));
end


function data = bitrotate(data,nBits)
% BITROTATE Bit-wise circular shift.
%    C = BITROTATE(A,K) returns the value of A circularly shifted by K bits
%    to the left. The left-most K bits that overflow as a result of the
%    shift are placed in the right-most K bits of C. A must be an unsigned
%    integer or a double (in which case A overflows and wraps around at the
%    implicit 53-bit limit for the floating point fraction of a double
%    precision number). A and K can be scalars or arrays of the same size.
%    Negative values of K are allowed and this corresponds to shifting to
%    the right and placing the right-most K overflowing bits in the
%    left-most K bits of C.
%
%    Example:
%       Circularly shift the bits of an unsigned 8-bit value to the left
%       until they return to their original positions. Display the bit
%       pattern for each shift.
%
%       a = uint8(63);
%       for i = 0:8
%          disp(dec2bin(bitrotate(a,i),8));
%       end
%
%    See also bitshift, bitand, bitor, bitxor, bitcmp, bitset, bitget.
% Author: Ken Eaton
% Version: MATLAB R2010b
% Last modified: 3/1/11
%--------------------------------------------------------------------------
  % Check the number and size of the inputs:
  assert(nargin == 2,'bitrotate:notEnoughInputs',...
         'Not enough input arguments.');
  assert(isscalar(nBits) || isscalar(data) || ...
         isequal(size(nBits),size(data)),'bitrotate:sizeMismatch',...
         'Inputs must have same size.');
  % Check the nBits input:
  if isempty(nBits)
    data = data([]);
    return
  end
  assert(all(rem(nBits(:),1) == 0),'bitrotate:inputsMustBeIntegers',...
         'Inputs must be integers.');
  % Check the data input:
  if isempty(data) || ~any(data(:))
    return
  end
  if isa(data,'double')
    assert(all(rem(data(:),1) == 0),'bitrotate:inputsMustBeIntegers',...
           'Inputs must be integers.');
    signData = sign(data);
    dataBits = 53;
  else
    [isUInt,typeIndex] = ismember(class(data),...
                                  {'uint8','uint16','uint32','uint64'});
    if isUInt
      dataBits = 2^(2+typeIndex);
    else
      error('bitrotate:badInputType',...
            'Input arguments must be of type unsigned integer or double.');
    end
  end
  % Reduce the number of bit rotations to a minimum (no need to rotate the
  %   bits by more than dataBits):
  nBits = rem(nBits,dataBits);
  if all(nBits == 0)
    return
  end
  % Find the portion of the data that overflows when shifted:
  overflowBits = (2.^abs(nBits)-1).*2.^((nBits > 0).*(dataBits-nBits));
  if isa(data,'double') && any(signData(:) < 0)
    overflowData = signData.*bitand(abs(data),overflowBits);
  else
    overflowData = bitand(data,overflowBits);
  end
  % Shift the data and add the overflow to the opposite end:
  data = bitshift(overflowData,nBits-sign(nBits)*dataBits)+...
         bitshift(data,nBits);
end