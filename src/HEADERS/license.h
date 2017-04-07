c+------------------------------------------- -*- fortran -*- ----------
c Distribution:  FEFF8.5L
c Copyright (c) [2013] University of Washington
c 
c This software was prepared in part with US Government Funding under
c DOE contract DE-FG03-97ER45623.

c Redistribution and use of this Distribution in source and binary
c formats, with or without modification is permitted, provided the 
c following conditions are met:
c 
c Redistributions must retain the above notices and the following list
c of conditions and disclaimer;
c 
c Modified versions carry the marking
c     "Based on or developed using Distribution: FEFF8.5L
c      Copyright (c) [2013] University of Washington"
c 
c Recipient acknowledges the right of the University of Washington to
c prepare uses of this Distribution and its modifications that may be
c substantially similar or functionally equivalent to
c Recipient-prepared modifications.
c
c Recipient and anyone obtaining access to the Distribution through
c recipient's actions accept all risk associated with possession and
c use of the Distribution.
c
c THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
c IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS TO THE
c DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c+------------------------------------------------------------------------------
c License is applicable for routines below, unless otherwise specified.
c
c+------------------------------------------------------------------------------
c Main contributors to feff8l, 2013-2017
c M. Newville,     newville AT cars DOT uchicago DOT edu
c B. Ravel,        bravel AT bnl DOT gov
c K. Jorissen,     kevinjorissen DOT pdx AT gmail DOT com
c J. J. Kas,       joshua DOT j DOT kas AT gmail DOT com
c
c+------------------------------------------------------------------------------
c Copyright Notice: FEFF8 is copyright protected software and users
c must obtain a license from the University of Washington Office of
c Technology Transfer for its use; see FEFF8 document for details.
c
c Main Authors of FEFF8: please contact us concerning any problems.
c A. L. Ankudinov
c B. Ravel,        bravel AT bnl DOT gov
c J. J. Rehr,      jjr AT uw DOT edu
c
c Citations: Please cite at least one of the following articles if 
c FEFF8 is used in published work: 
c    1) Main FEFF8 reference 
c       A.L. Ankudinov, B. Ravel, J.J. Rehr, and S.D. Conradson, 
c       Phys. Rev. B 58, 7565, (1998).
c       DOI: 10.1103/PhysRevB.58.7565
c    2) Multiple scattering theory
c       J.J. Rehr and R.C. Albers, Phys. Rev. B 41, 8139 (1990).
c       DOI: 10.1103/PhysRevB.41.8139
c+------------------------------------------------------------------------------
C
c Main Authors of FEFF5: please contact us concerning problems.
c A. L. Ankudinov
c S. I. Zabinsky
c J. J. Rehr,     jjr AT uw DOT edu
c R. C. Albers,   rca AT lanl DOT gov
C
c Citations: Please cite at least one of the following articles if 
c FEFF is used in published work: 
c    1) Multiple scattering
c       J.J. Rehr and R.C. Albers, Phys. Rev. B41, 8139 (1990).
c          DOI: 10.1103/RevModPhys.72.621
c       J.J. Rehr, S.I. Zabinsky and R.C. Albers, 
c          Phys. Rev. Let. 69, 3397 (1992).
c          DOI: 10.1103/PhysRevLett.69.3397
c    2) General reference
c       J.J. Rehr, J. Mustre de Leon, S.I. Zabinsky, and R.C. Albers,
c          J. Am. Chem. Soc. 113, 5135 (1991).
c          DOI: 10.1021/ja00014a001
c    3) Technical reference
c       J. Mustre de Leon, J.J. Rehr, S.I.  Zabinsky, and R.C. Albers,
c          Phys. Rev. B44, 4146 (1991).
c          DOI: 10.1103/PhysRevB.44.4146
c+------------------------------------------------------------------------------
