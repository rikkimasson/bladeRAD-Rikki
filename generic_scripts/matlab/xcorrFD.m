function [lags,ck,td] = xcorrFD(x,y)
%   FREQUENCY-DOMAIN CROSS-CORRELATION FUNCTION
%
%   Description: xcorrFD takes two discrete time signals as input and
%   calculates cross-correlation values and delay between two signals. The
%   computation is performed in the frequency domain. The results of
%   xcorrTD is validated against the MatLAB's xcorr function.
%
%   For cross-correlation in time domain see xcorrTD. 
%
%   Syntax:    
%         [ck,td] = xcorrFD(x,y)
%
%   Input: 
%         x = input signal 1 (must be a Nx1 vector)
%
%         y = input signal 2 (must be a Nx1 vector)
%
%   Output: 
%      lags = a vector of lags with a length of 2xN-1 (N = number of
%             data points in signal x or y)
%
%        ck = cross-correlation values [(2N-1)x1 vector] 
%
%        td = delay (i.e., number of lags) between two signals 
%
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
%   NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (ekalkan@usgs.gov)
%   $Revision: 1.1 $  $Date: 2017/11/27 14:03:00 $
a = y'; b = x';
len = length(a);
c = [ zeros(1,len-1) a ];
d = [ b zeros(1,len-1)  ];
% Compute FFTs
X1 = fft(c);
X2 = fft(d);
% Compute cross correlation 
X = X1.*conj(X2);
ck = ifft((X));
[~,i] = max(ck);
td = i - len;
lags = [-(len-1):len-1];