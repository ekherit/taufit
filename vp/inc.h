#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include"util.h"
#include"list.h"
#include"integral.h"
#include"init.h"
#include"R.h"
#include"lepton.h"
#include"hadron.h"
#include"const.h"

#define nr 8        // -- Число учтенных узких резонансов

double Mr[nr+1],Gr[nr+1],Beer[nr+1];

Tab TabR;
