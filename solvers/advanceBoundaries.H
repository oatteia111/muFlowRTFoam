// Advance once at the accepted interval start, never inside an equation/Picard loop.
tnext = GREAT;
if ((flowType > 0 && flowType < 4) || coupling == 1)
{
    advanceBoundary("hfix", presentTime, cellsHfix.size(), itHfix, tHfix, cellsHfixData);
    advanceBoundary("hwel", presentTime, cellsHwel.size(), itHwel, tHwel, cellsHwelData);
    advanceBoundary("hdrn", presentTime, cellsHdrn.size(), itHdrn, tHdrn, cellsHdrnData);
    advanceBoundary("hghb", presentTime, cellsHghb.size(), itHghb, tHghb, cellsHghbData);
    advanceBoundary("hriv", presentTime, cellsHriv.size(), itHriv, tHriv, cellsHrivData);
    advanceBoundary("hrch", presentTime, cellsHrch.size(), itHrch, tHrch, cellsHrchData);
}
if (flowType == 4)
{
    advanceBoundary("pfix", presentTime, cellsPfix.size(), itPfix, tPfix, cellsPfixData);
    advanceBoundary("pwel", presentTime, cellsPwel.size(), itPwel, tPwel, cellsPwelData);
    advanceBoundary("ffix", presentTime, cellsFfix.size(), itFfix, tFfix, cellsFfixData);
    advanceBoundary("fwel", presentTime, cellsFwel.size(), itFwel, tFwel, cellsFwelData);
}
if (activateTransport == 1 && activateReaction == 0 && activateEK == 0)
{
    advanceBoundary("chfix", presentTime, cellsChfix.size(), itChfix, tChfix, cellsChfixData);
    advanceBoundary("cfix", presentTime, cellsCfix.size(), itCfix, tCfix, cellsCfixData);
    advanceBoundary("cwel", presentTime, cellsCwel.size(), itCwel, tCwel, cellsCwelData);
    advanceBoundary("cghb", presentTime, cellsCghb.size(), itCghb, tCghb, cellsCghbData);
    advanceBoundary("criv", presentTime, cellsCriv.size(), itCriv, tCriv, cellsCrivData);
    advanceBoundary("crch", presentTime, cellsCrch.size(), itCrch, tCrch, cellsCrchData);
}
if (activateThermal == 1)
{
    advanceBoundary("thfix", presentTime, cellsThfix.size(), itThfix, tThfix, cellsThfixData);
    advanceBoundary("tfix", presentTime, cellsTfix.size(), itTfix, tTfix, cellsTfixData);
    advanceBoundary("twel", presentTime, cellsTwel.size(), itTwel, tTwel, cellsTwelData);
    advanceBoundary("tghb", presentTime, cellsTghb.size(), itTghb, tTghb, cellsTghbData);
    advanceBoundary("triv", presentTime, cellsTriv.size(), itTriv, tTriv, cellsTrivData);
    advanceBoundary("trch", presentTime, cellsTrch.size(), itTrch, tTrch, cellsTrchData);
}
if (activateReaction == 1)
{
    if (flowType < 4)
    {
        advanceBoundary("shfix", presentTime, cellsShfix.size(), itShfix, tShfix, cellsShfixData);
        advanceBoundary("sghb", presentTime, cellsSghb.size(), itSghb, tSghb, cellsSghbData);
        advanceBoundary("srch", presentTime, cellsSrch.size(), itSrch, tSrch, cellsSrchData);
    }
    else
    {
        advanceBoundary("spfix", presentTime, cellsSpfix.size(), itSpfix, tSpfix, cellsSpfixData);
    }
    advanceBoundary("sfix", presentTime, cellsSfix.size(), itSfix, tSfix, cellsSfixData);
    advanceBoundary("swel", presentTime, cellsSwel.size(), itSwel, tSwel, cellsSwelData);
    if (ph_gcomp > 0)
    {
        advanceBoundary("gfix", presentTime, cellsGfix.size(), itGfix, tGfix, cellsGfixData);
        advanceBoundary("gwel", presentTime, cellsGwel.size(), itGwel, tGwel, cellsGwelData);
    }
}
if (activateEK == 1)
{
    advanceBoundary("ewel", presentTime, cellsEwel0.size(), itEwel, tEwel, cellsEwelData);
}
