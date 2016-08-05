!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    a subroutine to increment the model and save it to model_data
!
!The MIT License (MIT)
!
!Copyright (c) 2016 Philip A. Browne
!
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!
!Email: p.browne@reading.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> subroutine to increment the model when the model is a subroutine
!! of empire. This is comms_v4 routine.
!! @param[out] model_states(:,particle)
subroutine model_as_subroutine_start(x,particle,tag)
  use model_as_subroutine_data
  use sizes
  implicit none
  real(kind=kind(1.0d0)), dimension(state_dim), intent(in) :: x !<
  !! the input state vector to be incremented
  integer, intent(in) :: particle           !< the particle number
  integer, intent(in) :: tag                !< the tag that controls
  !! whether to continue, finish or reset time
 
  !note that the output must be stored in model_states(:,particle)

  !below is an example of just the identity model. Please change this.
  model_states(:,particle) = x
  
end subroutine model_as_subroutine_start
