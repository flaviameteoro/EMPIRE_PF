!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    a subroutine to initialise and return the state from the model
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

!> subroutine to initialise and return the state from the model
subroutine model_as_subroutine_return(x,particle)
  use model_as_subroutine_data
  use sizes
  use comms, only : pfrank,nens,npfs
  implicit none
  real(kind=kind(1.0d0)), dimension(state_dim), intent(out) :: x
  integer, intent(in) :: particle
  integer :: i

  if(.not. initialised) then
     !allocate the space for the data
     first_ptcl = ceiling(real(pfrank)*real(nens)/real(npfs))+1
     final_ptcl = ceiling(real(pfrank+1)*real(nens)/real(npfs))
     allocate(model_states(state_dim,first_ptcl:final_ptcl))
     do i = first_ptcl,final_ptcl
        call model_as_subroutine_initialise(model_states(:,i),i)
     end do
  end if
  !now return the appropriate model state
  x = model_states(:,particle)
  
end subroutine model_as_subroutine_return
