        Integer n
        Character*80 ns

	Real*8 E,S,ReH,ImH,ReL,ImL

         n=IARGC()
         if(n.ne.1) then
          print *,'Usage: ./VP [Energy(GeV)]'
          return
         endif
         Call GetArg(1,ns)
         read(ns,*) E
         S=E**2

         Call vp_init  ! �������������� ��������� ������ � ����������
         Call vp_lepton(S,ReL,ImL) ! ������� ��������� ����� �����������
         Call vp_hadron(S,ReH,ImH) ! ������� �������� ����� �����������
         print *,E,ReL,ImL,ReH,ImH
         Call vp_done  ! ����������� ��������� ������ � ������, ������� ���

	End
