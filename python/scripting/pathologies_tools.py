from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring
import moby2
#from weather import APEX_weather
#import moby2.util.noninteractive_plots
import matplotlib, pylab
from matplotlib import pyplot as plt
from moby2.analysis.tod_ana import visual as v
pylab.ion()
from moby2.scripting import products
from moby2.analysis.tod_ana import pathologies
from moby2.util.database import TODList
import numpy, os, sys, time, ephem, pickle
np = numpy
pylab.ioff()




def recoverScanCuts(tod, params, full=False):
    """
    @brief  Script function that will find the cuts in a TOD and return them. If there is
            no saved version of the cuts, it will create them.

    Note that the pathologies object returned by this function is calibrated
    """
    cutParams = moby2.util.MobyDict.from_file(params.get("cutParams"))
    pathop = cutParams['pathologyParams']
    depot = moby2.util.Depot(params.get("depot"))
    name = tod.info.name
    if params.get("flatfield") is not None:
        cutParams['pathologyParams']['calibration']['flatfield'] = params.get("flatfield")


    # FIND STORED PATHOLOGIES
    pa = get_pathologies(tod, params)
    if pa is None:
        return None, None
    if tod.info.sample_index != pa.offsets[0] or tod.nsamps != pa.offsets[1]:
        fix_tod_length(tod, pa.offsets)

    # CUTS THAT DON'T NEED TO BE SAVED
    if cutParams.get('cutSlow',False):
        time_const = products.get_time_constants(cutParams.get('time_constants'), 
                                                 tod.info)
        f3dB = np.zeros(len(tod.det_uid))
        f3dB[time_const.det_uid] = time_const.tau
        pa.addCriterion('slow', f3dB > 1./(2*np.pi*cutParams.get('slowLevel',15)))
    if cutParams.get('cutBadRes',False):
        iv_cal = products.get_iv_calibration({"source":"data"}, tod.info)
        Rn_limits = params.get('fractionRn')
        Rn_pass = (Rn_limits[0] <= iv_cal.fraction_Rnormal) * \
                        (iv_cal.fraction_Rnormal < Rn_limits[1])
        pa.addCriterion('badRes', Rn_pass)

    # GET NEW PATHOLOGY SELECTIONS
    # Note that here we multiply all values by the flatfield or flatfield*responsivity
    pa.makeNewSelections(params=pathop)
    det_cuts = pa.makeCuts()

    # GENERATE CUTS OBJECT
    c_obj = moby2.TODCuts.for_tod(tod, assign = False)
    
    # APPLY ALL PARTIAL CUTS
    # MCE cuts
    c_obj.merge_tod_cuts(moby2.tod.get_mce_cuts(tod))
    # Glitch cuts
    cfp = depot.get_full_path(moby2.TODCuts, tag=params.get('tag_partial'), tod=tod)
    if os.path.exists( cfp ): 
        c_obj.merge_tod_cuts(depot.read_object( moby2.TODCuts, 
                                                tag=params.get("tag_partial"), 
                                                tod=tod))
    # Section cuts
    if pathop.get('getPartial',True):
        section_cuts = pa.makePartialCuts()
        c_obj.merge_tod_cuts(section_cuts)
    # Scan cuts
    if not(params.get('stare',True)): 
        az_cuts = pa.makeAzCuts()
        c_obj.merge_tod_cuts(az_cuts)
    
    # MERGE DETECTOR CUTS THAT DEPEND ON THE PARTIAL CUTS
    if cutParams.get('maxFraction') is not None:
        for det in tod.info.det_uid:
            frac = np.asarray(c_obj.cuts[det].get_mask(),dtype=int).sum()/float(c_obj.nsamps)
            if frac > cutParams['maxFraction']:
                det_cuts.set_always_cut(det)
                pa.liveSel[det] = False
    if cutParams.get_deep(("glitchParams","maxGlitch")) is not None:
        for det in tod.info.det_uid:
            if len(c_obj.cuts[det]) > cutParams["glitchParams"]["maxGlitch"]:
                det_cuts.set_always_cut(det)
                pa.liveSel[det] = False

    # MERGE DETECTOR CUTS
    c_obj.merge_tod_cuts(det_cuts)

    # ADD CALIRBATION CUTS (IF Calibration CANNOT BE COMPUTED, KILL WHOLE TOD)
    sel = pa.liveSel*pa.calData["respSel"]*(pa.crit["gainLive"]["values"] != 0)*pa.calData["stable"]
    if sel.sum() == 0:
        c_obj.set_always_cut(tod.info.det_uid)


    # GET DARK DETECTORS
    dd = pathologies.darkDets(patho = pa)

    # OUTPUT NEW CUTS
    if params.get('tag_out',None) is not None:
        assert params.get('tag_out') != params.get('tag_partial')
        depot.write_object(c_obj, tod = tod, force = True,
                           tag=params.get('tag_out')) 
        depot.write_object(det_cuts, tod = tod, force = True,
                           tag=params.get('tag_out')+"_det")
        depot.write_object( dd, tod = tod, force = True,
                            tag = params.get('tag_out'))
        if pathop.get('getPartial',True):
            depot.write_object(section_cuts, tod = tod, force = True,
                               tag=params.get('tag_out')+"_sec")
        print("cuts exported to %s"%params.get('tag_out'))

        if full:
            return c_obj, pa, resp, ff, calib
    return c_obj, pa


def get_pathologies(tod, params):
    cutParams = moby2.util.MobyDict.from_file(params.get("cutParams"))
    pathop = cutParams['pathologyParams']
    depot = moby2.util.Depot(params.get("depot"))
    fp = depot.get_full_path(pathologies.Pathologies, 
                             tag=params.get("tag_patho"), tod = tod,
                             structure = params.get("structure"))
    if os.path.exists(fp):
        pa = depot.read_object(pathologies.Pathologies, 
                               tag = params.get("tag_patho"),
                               tod = tod,
                               options = {"tod" : tod, "params": pathop},
                               structure = params.get("structure"))
        return pa
    else:
        print("Error: no pathology object found (as explicit argument or from depot)")
        return



# MAKE REPORT OF PATHOLOGIES AND CUTS RESULTS

class reportPathologies( object ):
    """
    @Brief Generate report files for pathology statistics
    """
    def __init__( self, params ):
        self.params = moby2.util.MobyDict.from_file(params)
        self.cutParams = moby2.util.MobyDict.from_file(self.params.get("cutParams"))
        p = self.params.get
        self.depot_file = os.path.abspath(os.path.join(p("outdir"),p("report")+".db"))

    def _initializeFiles( self, pa ):
        # ENTRIES FOR STATISTICS REPORT
        lkeys = pa.activeLiveKeys
        dkeys = pa.activeDarkKeys

        self.sel_keys = [ ("length", "float"),
                          ("liveDets","int"),
                          ("liveFrac","float"),
                          ("darkDets","int"),
                          ("glitches","int")]
        for k in lkeys:
            if "sel" in pa.crit[k]:
                self.sel_keys.append((k,"int"))
        for k in dkeys:
            if "sel" in pa.crit[k]:
                self.sel_keys.append((k,"int"))
        self.sel_keys.append(('gainCut','int'))
        self.sel_keys.append(('temperatureCut','int'))
        self.sel_keys.append(('Temp','float'))
        self.sel_keys.append(('DTemp','float'))

        self.stat_keys = []
        for k in lkeys:
            if "median" in pa.crit[k]:
                self.stat_keys.append((k,"float"))
        for k in dkeys:
            if "median" in pa.crit[k]:
                self.stat_keys.append((k,"float"))

        if self.params.get('newfile',False) or not(os.path.isfile(self.depot_file)):
            p = self.params.get
            rname = os.path.abspath(os.path.join(p("outdir"),p("report")))
            dname = os.path.dirname(rname)
            if not(os.path.isdir(dname)): os.makedirs(dname)

            hd1 = '# todName | '
            hd2 = '# str | '
            for k in self.sel_keys:
                hd1 += k[0] + ' | '
                hd2 += k[1] + ' | '
            for k in self.stat_keys:
                hd1 += '%s_m | %s_s | '%(k[0],k[0]) 
                hd2 += '%s | %s | '%(k[1],k[1])
            hd1 = hd1[:-3] + '\n'
            hd2 = hd2[:-3] + '\n'
         
            p = self.params.get
            cp = self.cutParams.get
            f = open(self.depot_file, 'w')
            f.write('# BEGIN HEADER\n')
            f.write('#\n# Date: %s\n' % time.asctime())
            f.write('#\n# tag: %s\n' % p('tag_patho'))
            f.write('# tag_partial: %s\n' % p('tag_partial'))
            f.write('#\n# Pathologies parameters:\n')
            f.write(printDictionary(cp('pathologyParams'), tabLevel = 1, 
                                    prefix = '#', verbose = False ))
            f.write('#\n# Glitches paramteres:\n')
            f.write(printDictionary(cp('glitchParams'), tabLevel = 1, 
                                    prefix = '#', verbose = False ))
            # f.write('#\n# Other cut parameters:\n')
            # f.write('#%28s = %d\n'%('cutBadRes',cp('cutBadRes')))
            # f.write('#%28s = [%12.4g ,%12.4g]\n' % ('fractionRn', 
            #         cp('fractionRn')[0], cp('fractionRn')[1]))
            # f.write('#%28s = %d\n'%('cutSlow',cp('cutSlow')))
            # f.write('#%28s = %12.4g Hz\n' % ('slowTimeC',cp('slowLevel')))
            f.write('# END HEADER\n#\n')
            f.write(hd1)
            f.write(hd2)
            f.close()

    def appendResult( self, obs ):
        """
        @Brief Add new entry to both results files.
        """
        p = self.params.get
        loadParams = {"filename":obs, "read_data":False}
        loadParams.update(p("params_loadtod",{}))
        tod = moby2.scripting.get_tod(loadParams)
        depot = moby2.util.Depot(p("depot"))
        pc_obj=depot.read_object(moby2.TODCuts, tag=p("tag_partial"), tod = tod)
        tod.info.sample_index = pc_obj.sample_offset
        glitches = len(tod.cuts.get_uncut())
        c_obj, pa = recoverScanCuts(tod, self.params)
        if "tag_cal" in self.params:
            _=recoverCalibration(tod,self.params, cuts=c_obj, pa=pa)
        self._initializeFiles(pa)
        # live = numpy.ones(len(pa.dets), dtype = 'bool')
        # live[1024:] = False

        # GET CUT STATISTICS
        darkDets = len(pa.dets[pa.darkSel])
        d = c_obj.get_uncut()
        liveDets = len(d)
        frac = 0.0
        for det in d:
            frac += numpy.asarray(c_obj.cuts[det].get_mask(),dtype=int).sum()/ \
                    float(c_obj.nsamps)/len(d)
        # GET TOD PARAMETERS
        length = (tod.ctime[-1] - tod.ctime[0])/60 # minutes
        Temp = pa.Temp
        DTemp = pa.dTemp

        # PRINT RESULTS
        f = open(self.depot_file, 'a')
        f.write("%s "%tod.info.name)
        f.write("%8.3g "%length)
        f.write("%4d "%liveDets)
        f.write("%8.3g "%frac)
        f.write("%4d "%darkDets)
        f.write("%4d "%glitches)
        for k in self.sel_keys[5:-4]:
            ndet = numpy.asarray(pa.crit[k[0]]["sel"]*~pa.origDark, dtype=int).sum()
            f.write('%4d '%ndet)
        f.write('%4d '%int(pa.gainCut))
        f.write('%4d '%int(pa.temperatureCut))
        if Temp is not None: f.write("%8.3g "%Temp)
        else: f.write("%8d "%(-1))
        if DTemp is not None: f.write("%8.3g "%DTemp)
        else: f.write("%8d "%(-1))
        for k in self.stat_keys:
            f.write('%8.3g %8.3g  '%(pa.crit[k[0]]['median'], pa.crit[k[0]]['sigma']))
        f.write('\n')
        f.close()


class pathoList( object ):
    """
    @brief Object that contains a list with results from a pathology run in the season
    """

    def __init__( self, filename, type = "pathoList" ):
        """
        @param type: "pathoList" or "TODList"
        """
        self.filename = filename
        if type == "pathoList": reader = readAscii
        elif type == "TODList": reader = readTODList
        else:
            psLib.trace("moby", 0, "ERROR: Unknown file type.")
            return -1

        if isinstance(filename, str): 
            self.data, self.header, self.keys, self.types = reader(filename)
        else: 
            self.data, self.header, self.keys, self.types = reader(filename[0])
            if len(filename) > 1:
                for i in range(1,len(filename)):
                    d, h, k, t = reader(filename[i])
                    for k in list(d.keys()): self.data[k] += d[k]

        self.selParams = {}
        self.keys.append("ctime")
        self.types.append("int")
        self.data['ctime'] = []
        for name in self.data['todName']:
            self.data['ctime'].append(int(name.split('/')[-1].split('.')[0]))
        self.ndata = len(self.data['todName'])
        for k in list(self.data.keys()):
            self.data[k] = numpy.array(self.data[k])


    def merge( self, pl2 ):
        """
        @Brief merge the self pathologies list with an external list keeping only unrepeated elements.
        """
        tod1 = self.data['todName']
        keys = list(self.data.keys())
        for i in range(pl2.ndata):
            if not(numpy.any(tod1 == pl2.data['todName'][i])):
                for k in keys:
                    self.data[k] = numpy.hstack([self.data[k],pl2.data[k][i]])
        s = numpy.argsort(self.data['todName'])
        for k in keys:
            self.data[k] = numpy.array(self.data[k])[s]
        self.ndata = len(s)


    def addPWV2results( self ):
        """
        """
        if "PWV" in self.keys: return
        w = get_pwv(self.data['ctime'])
        self.data["PWV"] = w
        self.keys.append("PWV")
        self.types.append("float")

    def addEphem( self ):
        """
        @brief Adds fields for the hour in the day of the observation and 
        hours from sunrise and to sunset.
        """
        act = ephem.Observer()
        act.lat = -(22+57./60+35./3600)*numpy.pi/180
        act.long = -(67+47./60+13./3600)*numpy.pi/180
        sun = ephem.Sun()
        self.data["hour"] = []
        self.data["hourAfterSunrise"] = []
        self.data["hourAfterSunset"] = []
        self.keys.extend(["hour","hourAfterSunrise","hourAfterSunset"])
        self.types.extend(["float", "float", "float"])
        for ct in self.data["ctime"]:
            act.date = ephem.Date(time.gmtime(ct)[:6])
            self.data["hour"].append((act.date-act.previous_antitransit(sun))*24)
            self.data["hourAfterSunrise"].append((act.date-act.previous_rising(sun))*24)
            self.data["hourAfterSunset"].append((act.date-act.previous_setting(sun))*24)
        self.data["hour"] = numpy.array(self.data["hour"])
        self.data["hourAfterSunrise"] = numpy.array(self.data["hourAfterSunrise"])
        self.data["hourAfterSunset"] = numpy.array(self.data["hourAfterSunset"])

    def addDB( self, year, array ):
        """
        @brief Adds extra information about each TOD obtained from the database
        """
        from moby2.instruments import actpol
        self.data["obs_type"] = []
        self.data["obs_detail"] = []
        self.data["obs_drift"] = []
        self.data["alt"] = []
        self.data["az_min"] = []
        self.data["az_max"] = []
        self.data["scan_speed"] = []
        self.keys.extend(["obs_type", "obs_detail", "obs_drift",
                          "alt, az_min", "az_max", "scan_speed"])
        self.types.extend(["str","str","str","float","float","float","float"])

        # Get database
        db = actpol.TODDatabase(config_file='/data/manifest_conf/manifest_%s.conf'%year)
        obs_type = ["planet","scan"]
        ids = db.select_tods(array=array, obs_type=obs_type)
        for id in self.data["todName"]:
            for r in ids:
                if r.basename == id:
                    self.data["obs_type"].append(r.obs_type)
                    self.data["obs_detail"].append(r.obs_detail)
                    self.data["obs_drift"].append(r.obs_drift)
                    self.data["alt"].append(min([r.mean_alt,61]))
                    self.data["az_min"].append(r.min_az)
                    self.data["az_max"].append(r.max_az)
                    if r.scan_speed is None or r.scan_speed < 0:
                        self.data["scan_speed"].append(1.5)
                    else:
                        self.data["scan_speed"].append(r.scan_speed)
        self.data["obs_type"] = np.array(self.data["obs_type"])
        self.data["obs_detail"] = np.array(self.data["obs_detail"])
        self.data["obs_drift"] = np.array(self.data["obs_drift"])
        self.data["alt"] = np.array(self.data["alt"])
        self.data["az_min"] = np.array(self.data["az_min"])
        self.data["az_max"] = np.array(self.data["az_max"])
        self.data["scan_speed"] = np.array(self.data["scan_speed"])




    def plotSeason( self, keywords, selection = None, factors = None, pwv = False, \
                    invSel = False, xlim = None, ylim = None, legLoc = 'upper right', \
                    xlabel = None, ylabel = None, grid = False, doubleAxis = False, \
                    title = None, filename = None, display = True):
        """
        """
        if len(keywords) != 2: doubleAxis = False
        if not(display): pylab.ioff()
 
        #date = []
        #for name in self.data['todName']:
        #    date.append(float(name.split('.')[0])/86400 + 2440587.5 - 4./24 - 1721424.5)
        #date = numpy.array(date)
        date = numpy.array(self.data['ctime'], dtype = float)/86400 + 2440587.5 - 4./24 - 1721424.5
 
        if factors is None: factors = numpy.ones(len(keywords))
        if selection is None: selection = numpy.ones(len(date), dtype = 'bool')
 
        if pwv:
            if doubleAxis: 
                pwv = False
                print("ERROR: cannot plot PWV in double axis mode. Omitting")
            else:
                ctime = numpy.array(self.data['ctime'], dtype = int)[selection]
                ct = (ctime - 54000)/86400
                night_ctime = []
                for c in range(ct.min(),ct.max()+1):
                    if numpy.any(ct == c):
                        night_ctime.append((c+1)*86400 + 3*3600)
                ct_ini = numpy.array(night_ctime).min()
                ct_end = numpy.array(night_ctime).max() + 86400
                ct_pwv = list(range(ct_ini, ct_end, 86400))
                PWV = get_pwv(ct_pwv)
                ct_date = []
                for c in ct_pwv:
                    ct_date.append(c/86400 + 2440587.5 - 1721424.5 - 3./24)
                #for i in range(len(PWV)):
                #    if PWV[i] == -1.0: PWV[i] = numpy.nan 

        mondays = pylab.WeekdayLocator(pylab.MO)
        months = pylab.MonthLocator()
        monthfmt = pylab.DateFormatter('%b %d')
        fo = matplotlib.font_manager.FontProperties(size=10)
 
        ax = pylab.subplot(111)
        lines = []
        leg = []
        if doubleAxis:
            lines += pylab.plot_date(date[selection], numpy.array(self.data[keywords[0]])[selection], \
                                     'b.', label = keywords[0])
            leg.append(keywords[0])
            if invSel:
                lines += pylab.plot_date(date[~selection], \
                                         numpy.array(self.data[keywords[0]])[~selection], \
                                         'b+', label = keywords[0]+" cut")
                leg.append(keywords[0])

            ax2 = pylab.twinx()
            ax2.yaxis.tick_right()
            lines += pylab.plot_date(date[selection], numpy.array(self.data[keywords[1]])[selection], \
                                     'g.', label = keywords[1])
            leg.append(keywords[1])
            if invSel:
                lines += pylab.plot_date(date[~selection], \
                                         numpy.array(self.data[keywords[1]])[~selection], \
                                         'g+', label = keywords[1]+" cut")
                leg.append(keywords[1])
        else:
            for i in range(len(keywords)):
                lines += pylab.plot_date(date[selection], \
                                         numpy.array(self.data[keywords[i]])[selection]*factors[i],'.')
                leg.append(keywords[i])
                if invSel:
                    lines += pylab.plot_date(date[~selection], \
                                   numpy.array(self.data[keywords[i]])[~selection]*factors[i], \
                                   '.', label = keywords[i]+" cut")
                    leg.append(keywords[i])
            if pwv:
                ax2 = pylab.twinx()
                ax2.yaxis.tick_right()
                pylab.plot(ct_date, PWV, 'k:')
                
 
        pylab.legend(lines, leg, numpoints = 1, prop = fo, loc = legLoc)
        #ax.xaxis.set_major_locator(mondays)
        ax.xaxis.set_major_formatter(monthfmt)
        #ax.xaxis.set_minor_locator(mondays)
        labels = ax.get_xticklabels()
        pylab.setp(labels, 'rotation', 30, fontsize=10)
        if grid: pylab.grid(True)
        ax.autoscale_view(tight = True, scaley = False)
 
        if xlim is not None: ax.set_xlim(xlim)
        if ylim is not None: ax.set_ylim(ylim)
        if xlabel is not None: ax.set_xlabel(xlabel)
        if doubleAxis:
            labels = ax2.get_xticklabels()
            pylab.setp(labels, 'rotation', 30, fontsize=10)
            ax2.set_frame_on(False)
            factor = factors[1]/factors[0]
            if xlim is not None: ax2.set_xlim(xlim)
            if ylim is not None: ax2.set_ylim([ylim[0]/factor,ylim[1]/factor])
            if ylabel is not None: 
                ax.set_ylabel(ylabel[0])
                ax2.set_ylabel(ylabel[1])
        if pwv:
            labels = ax2.get_xticklabels()
            pylab.setp(labels, 'rotation', 30, fontsize=10)
            ax2.set_frame_on(False)
            ax2.set_xlim(ax.get_xlim())
            ax2.set_ylim(0,4)
            ax2.set_ylabel('PWV [mm]')
        elif ylabel is not None: ax.set_ylabel(ylabel)
        if title is not None: pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if not(display): pylab.clf()
        else: return pylab.gca()


    def plotSeasonGrouped( self, keyword, selection = None, marker = None, color = None, \
                     invSel = False, includeCuts = False, xlim = None, ylim = None, \
                     xlabel = None, ylabel = None, pwv = False, weights = None, \
                     title = None, filename = None, display = True):
        """
        @brief   Plot statistic over season grouping daily results in a single point + error bars.
        @param  pwv   Add line indicating a smooth PWV contour.
        """
        if not(display): pylab.ioff()
        if selection is None: selection = numpy.ones(self.ndata, dtype = 'bool')
        elif invSel: selection = ~selection
        if weights is None: weights = numpy.ones(self.ndata, dtype = float)
        weights = weights[selection]
        data = numpy.array(self.data[keyword])[selection]
        ctime = numpy.array(self.data['ctime'], dtype = int)[selection]
        ct = (ctime - 54000)/86400
        if includeCuts:
            data_cut = numpy.array(self.data[keyword])[~selection]
            ctime_cut = numpy.array(self.data['ctime'], dtype = int)[~selection]
            date_cut = []
            for c in ctime_cut:
                date_cut.append(c/86400 + 2440587.5 - 1721424.5 - 3./24)
        
        night_mean = []
        night_std  = []
        night_ctime = []
        for c in range(ct.min(),ct.max()+1):
            day = (ct == c)
            if numpy.any(day):
                w = weights[day]/weights[day].sum()
                m = numpy.sum(data[day]*w)
                s = numpy.sqrt(numpy.sum(numpy.power(data[day]-m,2)*w))
                night_mean.append(m)
                night_std.append(s)
                night_ctime.append((c+1)*86400 + 3*3600)

        if pwv: 
            ct_ini = numpy.array(night_ctime).min()
            ct_end = numpy.array(night_ctime).max() + 86400
            ct_pwv = list(range(ct_ini, ct_end, 86400))
            PWV = get_pwv(ct_pwv)
            ct_date = []
            for c in ct_pwv:
                ct_date.append(c/86400 + 2440587.5 - 1721424.5 - 3./24)
            #for i in range(len(PWV)):
            #    if PWV[i] == -1.0: PWV[i] = numpy.nan 
                

        date = []
        for c in night_ctime:
            date.append(c/86400 + 2440587.5 - 1721424.5 - 3./24)


        ax = pylab.subplot(111)
        if marker is None: marker = 'o'
        if color is None: color = 'k'

        if pwv:
            pylab.errorbar(date, night_mean, yerr = night_std,
                        color = color, ls = ' ', marker = marker, mfc = color)                                   
            if includeCuts: pylab.plot(date_cut, data_cut, 'kx')
            ax2 = pylab.twinx()
            ax2.yaxis.tick_right()
            pylab.plot(ct_date, PWV, '%s--'%color)
        else:
            pylab.errorbar(date, night_mean, yerr = night_std,
                        color = color, ls = ' ', marker = marker, mfc = 'w')                                   
            if includeCuts: pylab.plot(date_cut, data_cut, 'kx')

        monthfmt = pylab.DateFormatter('%b %d')
        ax.xaxis.set_major_formatter(monthfmt)
        labels = ax.get_xticklabels()
        pylab.setp(labels, 'rotation', 30, fontsize=14)
        ax.autoscale_view(tight = True, scaley = False)
        if xlim is not None: ax.set_xlim(xlim)
        if ylim is not None: ax.set_ylim(ylim)
        if xlabel is not None: ax.set_xlabel(xlabel)
        else: pylab.xlabel('Date')
        if pwv:
            labels = ax2.get_xticklabels()
            pylab.setp(labels, 'rotation', 30, fontsize=14)
            ax2.set_frame_on(False)
            ax2.set_xlim(ax.get_xlim())
            ax2.set_ylim(0,4)
            ax2.set_ylabel('PWV [mm]')
        if ylabel is not None: ax.set_ylabel(ylabel)
        else: ax.set_ylabel(keyword)
        if title is not None: pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if not(display): pylab.clf()
        else: return pylab.gca()
        

    def histogram(self, data, selection = None, bins = None, range = None, \
                  title = None, xlabel = None, ylabel = None, \
                  filename = None, display = True, normed = False, color = "w", alpha = 1.):
        """
        @brief   Make histogram for the value for various statistics.
        @param   data      data vector (ndet long) to plot.

        @param  selection bool array with selection of detectors to include in the plot.
        @param  bins      number of bins to use.
        @param  title     string with title to use in figure.
        @param  filename  string with name of the file where to save the figure.
        @param  display   display or not the image.
        """
        if not(display): pylab.ioff()
        if selection is None:
            selection = numpy.ones(len(data), dtype = 'bool')
        dat = numpy.array(data)[selection].flatten()

        if range is not None:
            dat = dat[(dat >= range[0])*(dat <= range[1])]

        if bins is None:
            bins = int(len(dat)/100)

        a = pylab.hist(dat, bins = bins, fc = color, normed = normed, alpha = alpha)
        if xlabel is not None: pylab.xlabel(xlabel)
        if ylabel is not None: pylab.ylabel(ylabel)
        if title is not None: pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if display: return pylab.gca()
        else: pylab.clf()


    def getSignificance( self, key, increasing = True, selection = None, filename = None, steps = 10,
                         display = True):
        """
        @brief Obtain plot and table with percentage of usable data lost by a given file selection criterion
        @params   key        criterion to analyse
        @params   increasing wether to plot for ingreasing (gt) or decreasing (lt) values of key
        @params   selection  pre-selection used before analyzing one criterion
        @params   filename   file to output (produces a .png plot and a .txt table)
        @params   steps      how many steps to use
        @params   show       wether to show the plot
        """
        if not(display): pylab.ioff()
        assert numpy.any(numpy.array(list(self.data.keys())) == key)
        if ~numpy.any(numpy.array(list(self.data.keys())) == 'Effective'):
            print("ERROR: no 'Effective' field found. You can generate it as 'liveDets'*'fracTimeLive'")
            return
        if selection is None: selection = numpy.ones(self.ndata, dtype = bool)
        eff = self.data['Effective'][selection]
        total = eff.sum()
        data = self.data[key][selection]
        exp = numpy.floor(numpy.log10(data.max()))
        max = numpy.ceil(data.max() / numpy.power(10,exp)*2)/2*numpy.power(10,exp)
        points = numpy.arange(steps+1)*max/steps
        result = []
        for r in points:
            if increasing:
                result.append(eff[data < r].sum()/total * 100)
            else:
                result.append(eff[data > r].sum()/total * 100)

        if filename:
            if len(filename.split('.')) > 1:
                root = filename.split('.')[:-1]
            else:
                root = filename
            f = open('%s.txt'%root, 'w')
            f.write('# %s | Cut %s\n'%(key,'%'))
            for i in range(len(points)):
                f.write('%8.1f %8.2f\n'%(points[i], result[i]))
            f.close()
        
        pylab.plot(points, result, '-o')
        pylab.xlabel('Cut Level (%s)'%key)
        pylab.ylabel('Percentage of data cut')
        if filename: pylab.savefig('%s.png'%root)
        if display: return pylab.gca()
        else: pylab.clf()


    def reportStatistics( self, selection = None, filename = None, verbose=False):
        """
        @brief Obtain mean and standard deviation data in pathoList
        @param filename  output the results to a text file
        """
        if selection is None: selection = numpy.ones(self.ndata, dtype = bool)
        if filename:
            f = open(filename, 'w')
            f.write('# PathoList statistics (mean and standard deviation)')
        l = '# %-20s %+12s %+12s %+12s'%('Key', 'Mean', 'Median', 'Std')
        if verbose: print(l)
        if filename: f.write('\n%s\n'%l)
        for k in self.keys:
            if k != 'todName':
                s = ~numpy.isnan(self.data[k])
                l = '%-20s: %12.4g %12.4g %12.4g'%(k,
                                self.data[k][selection*s].mean(),
                                np.median(self.data[k][selection*s]),
                                self.data[k][selection*s].std())
                if verbose: print(l)
                if filename:
                    f.write('%s\n'%l)
        l = 'Number of TODs: %d'%len(numpy.where(selection)[0])
        if "length" in self.keys:
            obstime = self.data["length"][selection].sum()/60.
            l += "\nObserving time: %12.1f hours"%obstime
            if ("liveDets" in self.keys) and ("liveFrac" in self.keys):
                dettime = (self.data["length"][selection]*\
                          self.data["liveDets"][selection]*\
                          (1. - self.data["liveFrac"][selection])).sum()/60.
                effdets = dettime/obstime
                l += "\n1 detector time: %12.1f hours"%dettime
                l += "\nEffective detectors: %12.1f"%effdets
        if verbose: print(l)
        if filename:
            f.write('\n%s\n'%l)




    def selectFromTODList( self, filename, verbose = False ):
        """
        @brief Obtain a selection for the pathoList matching the TODs from
               a file
        @param   filename   TOD list filename
        @param   verbose    print extra information
        @return  selection  Bool array with selected TODs
        """
        sel = numpy.zeros(self.ndata, dtype = bool)
        todNames = numpy.array(self.data['todName'])
        f = open(filename, 'r')
        for l in f:
            if l[0] != '#':
                name = l.split('\n')[0].split('#')[0].split()[0].split('/')[-1].strip()
                sel[todNames == name] = True
                if verbose:
                    if ~numpy.any(todNames == name):
                        print('TOD %s not found')
        return sel

    
    def applySelection(self, selection):
        """
        @brief leave only the selected entries in object
        """
        self.ndata = selection.sum()
        for k in self.keys:
            self.data[k] = self.data[k][selection]


    def selectGoodTODs(self, params = None, preselection = None):
        """
        @brief  Generates a list of TODs that satisfy a set of requirements
        @params params      Dictionary with keys and conditions to search for ('lt': <, 'gt': >)
        @return selection   Bool array with selected TODs
        """
        if params is not None:
            self.selParams.update(params)

        selection = numpy.ones(len(self.data['todName']), dtype = 'bool')
        if preselection is not None: selection *= preselection

        keys = list(self.selParams.keys())
        for k in keys:
            if not isinstance(self.selParams[k],dict): continue
            if k in self.keys:
                if "lt" in list(self.selParams[k].keys()):
                    selection *= (self.data[k] < self.selParams[k]['lt'])
                if "gt" in list(self.selParams[k].keys()):
                    selection *= (self.data[k] > self.selParams[k]['gt'])
            else:
                print('WARNING: Selection key %s not found' % k)

        return selection

    def matchTODs(self, other):
        """
        @brief  Find TOD entries that match this list with another list.
        @param  other   Another pathoList object or TOD name list
        """
        if isinstance(other, pathoList): tods = other.data["todName"]
        else: tods = other
        mine = []
        theirs = []
        for o in range(len(tods)):
            i = numpy.where(self.data["todName"] == tods[o])[0]
            if len(i) > 0:
                mine.append(i[0])
                theirs.append(o)
        return mine, theirs

    def removeIndex(self, index):
        """
        @brief Removes a given data element
        @params     index   element to remove
        """
        assert index < self.ndata
        for k in self.keys:
            self.data[k] = np.hstack([self.data[k][:index],self.data[k][index+1:]])
        self.ndata -= 1


    def removeDuplicates(self, verbose = False):
        """
        @brief removes all duplicate elements from data
        """
        i = 0
        c = 0
        while i < len(self.data["todName"]):
            inds = np.where(self.data["todName"]==self.data["todName"][i])[0]
            if len(inds) > 1:
                if verbose:
                    for j in inds:
                        print("%s %d %d"%(self.data["todName"][j],
                                          self.data["liveDets"][j],
                                          self.data["darkDets"][j]))
                for j in inds[1:]:
                    self.removeIndex(j)
                c += 1
            i += 1
        print("Removed %d duplicates"%c)


    def outputTODList(self, filename, selection = None):
        """
        @brief  Outputs a TOD list for the specified selection
        @params   filename   output filename
        @params   selection  selected data points
        @params   addHeader  wether to add the header to the output
        """
        fb = moby2.scripting.get_filebase()

        if selection is None: 
            selection = numpy.ones(self.ndata, dtype = 'bool')

        header = '# BEGIN HEADER\n'
        header += '#\n# TOD selection parameters:\n'
        for k in list(self.selParams.keys()):
            if not isinstance(self.selParams[k],dict): continue
            header += '#    %20s:'%k
            for kk in list(self.selParams[k].keys()):
                    header += '    %30s = %12.4g\n' % (kk,float(self.selParams[k][kk]))
        header += '#\n# END HEADER\n'
        selTODs = self.data["todName"][selection].tolist()
        tl = TODList([fb.filename_from_name(name)[0] for name in selTODs])
        tl.to_file(filename, header = header)

    def outputSubset(self, filename, keys = None, selection = None, addHeader = True):
        """
        @brief  Outputs a pathoList with the specified fields and selection
        @params   keys       list of fields to include
        @params   filename   output filename
        @params   selection  selected data points
        @params   addHeader  wether to add the header to the output
        """
        if keys is None: keys = self.keys
        else:
            ind = (numpy.array(keys) == 'todName')
            if not(numpy.any(ind)):
                keys = ['todName'] + keys
            else:
                keys = ['todName'] + list(numpy.array(keys)[~ind])

        if selection is None: 
            selection = numpy.ones(len(self.data['todName']), dtype = 'bool')
            addHeader = False
        n = len(numpy.array(self.data['todName'])[selection])

        hd1 = '# '
        for k in keys:
            hd1 += k + ' | '
        hd1 = hd1[:-3] + '\n'

        hd2 = '# '
        for k in keys:
            ty = type(self.data[k][0])
            if ty == str or ty == numpy.string_: 
                hd2 += 'str'
            elif ty == float or ty == numpy.float64 or ty == numpy.float32:
                hd2 += 'float'
            elif ty == int or ty == numpy.int64 or ty == numpy.int32:
                hd2 += 'int'
            hd2 += ' | '
        hd2 = hd2[:-3] + '\n'
         
        f = open(filename, 'w')
        f.write('# BEGIN HEADER\n')
        for h in self.header: f.write('#%s\n' % h)
        if addHeader:
            f.write('#\n# TOD selection parameters:\n')
            for k in list(self.selParams.keys()):
                f.write('#    %20s:\n' % k)
                for kk in list(self.selParams[k].keys()):
                    f.write('#    %30s = %8.3g\n' % (kk,float(self.selParams[k][kk])))
        f.write('#\n')
        f.write('# END HEADER\n')
        f.write(hd1)
        f.write(hd2)
        for i in range(n):
            for k in keys:
                ty = type(self.data[k][0])
                if ty == str or ty == numpy.string_:
                    f.write(numpy.array(self.data[k])[selection][i] + ' ')
                elif ty == float or ty == numpy.float64 or ty == numpy.float32:
                    f.write('%8.3g' % numpy.array(self.data[k])[selection][i] + ' ')
                elif ty == int or ty == numpy.int64 or ty == numpy.int32:
                    f.write('%4d' % numpy.array(self.data[k])[selection][i] + ' ')
            f.write('\n')
        f.close()


    
def readAscii( filename ):
    header = []
    f = open(filename)
    names = f.readline()
    names = names.split("#")[1]
    names = names.split("\n")[0]
    if names.strip() == 'BEGIN HEADER':
        while names.strip() != 'END HEADER':
            names = f.readline()
            names = names.split("#")[1]
            names = names.split("\n")[0]
            header.append(names)
        del header[-1]
        names = f.readline()
        while len(names.split()) == 1: names = f.readline()
        names = names.split("#")[1]
        names = names.split("\n")[0]
    names = names.split("|")
    ft = f.readline()
    ft = ft.split("#")[1]
    ft = ft.split("\n")[0]
    ft = ft.split("|")
    assert len(names) == len(ft)
    frmt = {}
    data = {}
    names_output = []
    ft_output = []
    for i in range(len(names)):
        frmt[names[i].strip()] = {'type':ft[i].strip(), 'column':i}
        data[names[i].strip()] = []
        names_output.append(names[i].strip())
        ft_output.append(ft[i].strip())
    for l in f:
        ll = l.split()
        for k in list(frmt.keys()):
            data[k].append(eval(frmt[k]["type"])(ll[frmt[k]["column"]]))
    f.close()
    return data, header, names_output, ft_output

def getdtype(st):
    s1 = st.split("-")
    if (len(s1) == 1) or (s1[0] == ""):
        if s1[-1].isdigit(): return int
        s2 = s1[-1].split(".")
        if len(s2) == 2:
            if s2[0].isdigit() and s2[1].isdigit(): return float
    return str
                                                

def readTODList( filename ):
    header = []
    names_output = ["todName"]
    ft_output = ["str"]
    f = open(filename)
    obsname = []
    for l in f:
        if l[0] != "#":
            obsname.append(l.split("\n")[0].split()[0].split("/")[-1])
    data = {"todName": numpy.array(obsname)}
    return data, header, names_output, ft_output

def printDictionary( dictio, tabLevel = 0, tabSize = 4, prefix = '', verbose = False ):
    """
    @brief Visualize dictionary contents
    """
    output = ''
    for k in list(dictio.keys()):
        if isinstance(dictio[k], dict):
            output += "%s%s%-20s:\n"%(prefix,tabLevel*tabSize*' ', k)
            output += printDictionary(dictio[k], prefix = '#', tabLevel = tabLevel+1)
        elif isinstance(dictio[k], float):
            output += "%s%s%-20s: %12.3g\n"%(prefix,tabLevel*tabSize*' ', k, dictio[k])
        elif isinstance(dictio[k], int):
            output += "%s%s%-20s: %12d\n"%(prefix,tabLevel*tabSize*' ', k, dictio[k])
        elif isinstance(dictio[k], list):
            l = '['
            for i in range(len(dictio[k])):
                if i > 0: l += ' '
                if isinstance(dictio[k][i], int): l += '%d'%dictio[k][i]
                elif isinstance(dictio[k][i], float): l += '%g'%dictio[k][i]
                else: l += ' %s'%dictio[k][i]
                if i < len(dictio[k])-1: l += ','
            l += ']'
            output += "%s%s%-20s: %12s\n"%(prefix,tabLevel*tabSize*' ', k, l)
        else:
            output += "%s%s%-20s: %12s\n"%(prefix,tabLevel*tabSize*' ', k, dictio[k])
    if verbose: print(output)       
    return output

def get_pwv(ctimes, saturate = True, time_diff = 3600.):
    """                                                                        
    @brief Obtain PWV from APEX log files stored in the specified directory.   
    """
    apex_pwv = moby2.aux_data.apex.WeatherChannel('radiometer')
    w, tdiff = apex_pwv.get_nearest(ctimes)
    w[np.abs(tdiff)>time_diff] = -1.
    if saturate: w[w>10] = 10
    return w

def fix_tod_length(tod, offsets):
    tod.info.sample_index = offsets[0]
    tod.nsamps = offsets[1]
    tod.ctime = tod.ctime[offsets[0]:offsets[0]+offsets[1]]
    tod.az = tod.az[offsets[0]:offsets[0]+offsets[1]]
    tod.alt = tod.alt[offsets[0]:offsets[0]+offsets[1]]
    tod.data_mask = tod.data_mask[offsets[0]:offsets[0]+offsets[1]]
    if tod.data is not None:
        tod.data = tod.data[:,offsets[0]:offsets[0]+offsets[1]].copy()
    tod.enc_flags = tod.enc_flags[offsets[0]:offsets[0]+offsets[1]]


###############################################################################
####################        CALIBRATION TOOLS        ##########################
###############################################################################


def recoverCalibration(tod, params, pa=None, cuts=None, **kwargs):
    cutParams = moby2.util.MobyDict.from_file(params.get("cutParams"))
    pathop = cutParams['pathologyParams']
    depot = moby2.util.Depot(params.get("depot"))
    name = tod.info.name
    if params.get("flatfield") is not None:
        pathop['calibration']['flatfield'] = params.get("flatfield")

    if pa is None:
        pa = get_pathologies(tod, params)
        if pa is None: 
            return
    if tod.info.sample_index != pa.offsets[0] or tod.nsamps != pa.offsets[1]:
        fix_tod_length(tod, pa.offsets)

    if cuts is None:
        fp = depot.get_full_path(moby2.TODCuts,
                                 tag=params.get("tag_out"), tod = tod,
                                 structure = params.get("structure"))
        if os.path.exists(fp):
            cuts = depot.read_object(moby2.TODCuts, 
                                      tag = params.get("tag_out"),
                                      tod = tod,
                                      structure = params.get("structure"))
        else:
            print("Error: no cuts object found (as explicit argument or from depot)")
            return

    # GET CALIBRATION WITH ATM CORRECTION
    resp, ff, ffRMS, respSel, ffSel, stable = pa.getpWCalibration()

    # Force the pathologies object to be calibrated
    # This implies that the gains have already been multipied by the flatfield
    # meaning that they correspond to the gain excess after flatfielding
    if not(pa.calibrated): pa.calibrateValues()
    gains = pa.crit["gainLive"]["values"]

    # Select detectors that passed the cuts and with valid responsivities
    # and which are stable. Also make sure they are nonzero
    liveSel = cuts.get_mask()
    sel = liveSel*respSel*stable#*(gains != 0)
    if sel.sum() == 0:
        print("Unable to calibrate")
        return 0

    # Set the level
    if "level_type" not in pathop["calibration"]:
        pathop["calibration"]["level_type"] = "good"

    calib = np.zeros(pa.dets.size)
    mask = np.zeros(pa.dets.size)
    level = np.zeros(pa.dets.size)
    freqs = np.unique(tod.info.array_data.get('nom_freq',["single"]))

    for freq in freqs:
        if freq == "single":
            sf = np.ones(pa.dets.size,dtype=bool)
        else:
            sf = tod.info.array_data['nom_freq'] == freq
        if freq in pathop["calibration"].get("noATMfreq",[]):
            calib[sf] = resp[sf]*ff[sf]
        else:
            if sel[sf].sum() > 0:
                calib[sf], mask[sf], level[sf] = compute_calibration(
                                             gains[sf],ff[sf],resp[sf],sel[sf],
                                             **pathop['calibration'])

    calib *= tod.info.array_data["optical_sign"]
    # Store results
    if "tag_cal" in params:
        s = pa.liveSel 
        calObj = moby2.Calibration(det_uid = pa.dets[s])
        calObj.set_property(["cal", "calRMS"], [calib[s], ffRMS[s]])
        depot.write_object( calObj,
                            tag=params.get("tag_cal"),
                            tod = tod,
                            force = True)

    return calObj, mask, level


def compute_calibration(gains, ff, resp, sel,
                        level_type="good",**kwargs):

    if "min_sel" not in kwargs: kwargs["min_sel"] = 100
    if "min_stable" not in kwargs: kwargs["min_stable"] = 40

    if level_type is "good":
        weights = None
    elif level_type is "stable":
        weights = resp*ff
    else:
        raise ValueError("Error: Unknown level_type (should be either good or stable)")

    level, mask = pathologies.normalizeGains(gains, sel=sel, 
                                             weights = weights, **kwargs)

    nonzero = (gains != 0)
    cal = resp * ff
    cal[nonzero] /= gains[nonzero]
    return cal, mask, level


def cleanStable(sel, atm, delta = 0.2):
    """
    Selects only detectors with gains a delta fraction away from the median atm factor
    """
    m = np.median(atm[sel])
    n_atm = atm / m
    new_sel = sel*(n_atm>1.-delta)*(n_atm<1.+delta)
    return new_sel

def get_slope(x,y,xlim=[0,2.5],perc=90,plot=False,**kwargs):
    xsel = (x>xlim[0])*(x<xlim[1])
    ys = np.sort(y[xsel])
    err_sel = (y>ys[int(len(ys)*(1-perc/100.))])*(y<ys[int(len(ys)*(perc/100.))])
    for it in range(5):
        slope,y0 = np.polyfit(x[xsel*err_sel],y[xsel*err_sel],1)
        model = x*slope+y0
        err = y-model
        err_sel = (abs(err) < err.std())
#        err_sel = (abs(err) < err[err_sel].std())
    err5 = err[ abs(err)<err[err_sel].std()*5 ].std()
    err_near = err[err_sel].std()
    if plot:
        plot_noise_vs_sky(x,y,slope,y0,err5,err_near,**kwargs)
    return [slope, y0, err5, err_near]


def plot_noise_vs_sky(sky,val,slope,y0,err_all,err_near,
                      fign = 1,
                      xlim=[0,8],y_lim=None,
                      ylabel=None, title=None,
                      filename=None,show=False):
    pylab.ioff()
    x_lim = np.array([0,8])
    fit_ends = x_lim*slope+y0
    if y_lim is None: y_lim = [0,fit_ends[1]*2]
    print(fign)
    fig=pylab.figure(fign)
    pylab.plot(sky,val,".")
    pylab.plot(x_lim,fit_ends,"r")
    ax = fig.axes[0]
    pylab.xlim(x_lim)
    pylab.ylim(y_lim)
    pylab.text(0.5,0.9,"y = %8gx+%8g"%(slope,y0),transform=ax.transAxes)
    pylab.text(0.5,0.85,"err_rms = %8g"%err_all,transform=ax.transAxes)
    pylab.text(0.5,0.8,"err_rms (near) = %8g"%err_near,transform=ax.transAxes)
    pylab.xlabel("Sky [PWV/sin(alt)]")
    if ylabel is not None: pylab.ylabel(ylabel)
    if title is not None: pylab.title(title)
    if filename is not None: pylab.savefig(filename)
    if show:
        pylab.ion()
        pylab.show()
    else:
        pylab.close(fig)
        pylab.ion()

def read_pickle(filename):
    f = open(filename)
    p = pickle.Unpickler(f)
    data = p.load()
    f.close()
    return data

def get_metadata(names, array, year):
    from moby2.instruments import actpol
    alt = []
    ctime = []
    az_min = []
    az_max = []
    scan_speed = []
    obs_type = []
    obs_detail = []
    obs_drift = []
    m = {}
    obs_types = ["planet","scan"]
    db = actpol.TODDatabase(config_file='/data/manifest_conf/manifest_%s.conf'%year)
    ids = db.select_tods(array=array, obs_type=obs_types)
    for id in names:
        for r in ids:
            if r.basename == id:
                ctime.append(r.ctime)
                obs_type.append(r.obs_type)
                obs_detail.append(r.obs_detail)
                obs_drift.append(r.obs_drift)
                alt.append(min([r.mean_alt,61]))
                az_min.append(r.min_az)
                az_max.append(r.max_az)
                if r.scan_speed is None or r.scan_speed < 0:
                    scan_speed.append(1.5)
                else:
                    scan_speed.append(r.scan_speed)
    m["obs_type"] = np.array(obs_type)
    m["obs_detail"] = np.array(obs_detail)
    m["obs_drift"] = np.array(obs_drift)
    m["alt"] = np.array(alt)*np.pi/180.
    m["az_min"] = np.array(az_min)
    m["az_max"] = np.array(az_max)
    m["scan_speed"] = np.array(scan_speed)
    m["ctime"] = np.array(ctime)
    m["pwv"] = get_pwv(ctime)
    return m

    
    
def get_stable(data,meta,ff,N_min=500,outdir=None,it=0,
               respDisp=1.2e-5,calDisp=5e-5,atmDisp=0.1, minGain=0.2):

    # Get responsivities and calibrations
    shape = data["gainLive"].shape
    stable = ff.get_property("stable",det_uid=np.arange(shape[0]))[1]
    resp = np.array(data["resp"])
    nonzero = data["gainLive"] >= minGain
    sels = np.array(data["sel"],bool)*nonzero
    atm = np.ones(shape)
    atm[nonzero] = 1./data["gainLive"][nonzero]
    norms,goodTODs = normalizeAtms(atm,sels,stable)
    sels *= numpy.repeat([goodTODs],shape[0],axis=0)
    
    # Obtain calibration
    # ffc = ff.get_property("cal",det_uid=np.arange(shape[0]))[1]
    # ffc = numpy.repeat([ffc],shape[1],axis=0).T
    cal = np.ones(shape)
    cal[nonzero] = atm[nonzero]*resp[nonzero]#*ffc[nonzero]

    # Find detectors mostly uncut
    N = sels.sum(axis=1)
    G = np.where(N>N_min)[0]
    _,ffc = ff.get_property("cal",det_uid=G)

    # Get slopes
    slope_BSresp=[]
    slope_cal=[]
    slope_atm=[]

    id=""
    for i,fc in zip(G,ffc):
        s = sels[i]
        sky = meta["pwv"][s]/np.sin(meta["alt"][s])
        # RMS BSresp
        val = data["rmsLive"][i][s]*resp[i][s]*fc
        slope_BSresp.append(get_slope(sky,val))
        # RMS cal
        val = data["rmsLive"][i][s]*cal[i][s]
        slope_cal.append(get_slope(sky,val))
        # Atmosphere slope
        val = atm[i][s]
        slope_atm.append(get_slope(sky,val))


    slope_BSresp = np.array(slope_BSresp)
    slope_cal = np.array(slope_cal)
    slope_atm = np.array(slope_atm)


    # Find stable detectors
    stable = np.zeros(shape[0],dtype=bool)
    slope_cal_m, slope_cal_s = medsig(slope_cal[:,0])
    slope_BSresp_m, slope_BSresp_s = findMode(slope_BSresp[:,0],window=1.2e-5)
    sel = (slope_atm[:,0]>-atmDisp)*\
          (slope_atm[:,0]< atmDisp)*\
          (slope_cal[:,0]>slope_cal_m-2*slope_cal_s)*\
          (slope_cal[:,0]<slope_cal_m+2*slope_cal_s)*\
          (slope_BSresp[:,0]>slope_BSresp_m-respDisp/2)*\
          (slope_BSresp[:,0]<slope_BSresp_m+respDisp/2)*\
          (slope_cal[:,3]<calDisp)
    stable[G[sel]] = True
    
    # Make histogram of the slopes
    if outdir is not None:
        h1 = pylab.hist(slope_cal[:,0]*1e6,bins=250,alpha=0.5,range=[-50,200])
        h2 = pylab.hist(slope_BSresp[:,0]*1e6,bins=250,alpha=0.5,
                        range=[-50,200])
        h3 = pylab.hist(slope_BSresp[:,0][stable[G]]*1e6,bins=250,
                        alpha=0.5,range=[-50,200])
        pylab.legend(["CAL slope","BSresp slope","BSresp slope (stable)"])
        pylab.xlabel("Slope [aW/mm]")
        pylab.ylabel("# DETs")
        pylab.savefig("%s/hist_slope_comp_it%d.png"%(outdir,it))
        pylab.close()
        # Atmosphere slope histogram
        h1 = pylab.hist(slope_atm[:,0],bins=100,range=[-0.2,0.2])
        pylab.xlabel("Slope [1/mm]")
        pylab.ylabel("# DETs")
        pylab.savefig("%s/hist_slope_atm_it%d.png"%(outdir,it))
        pylab.close()

    # Find new flatfield
    norms,goodTODs = normalizeAtms(atm,sels,stable)
    sels *= numpy.repeat([goodTODs],shape[0],axis=0)
    N = sels.sum(axis=1)
    G = np.where(N>100)[0]
    ff_new = moby2.detectors.FlatField()
    ffcn=[];ffcne=[]
    for i in G:
        m,s=medsig(atm[i][sels[i]])
        ffcn.append(m)
        ffcne.append(s)
    ff_new.cal = np.array(ffcn)
    ff_new.calRMS = np.array(ffcne)
    ff_new.det_uid = G
    ff_new.stable = stable[G]

    return stable, ff_new
    

def normalizeAtms(atms,sels,stable):
    norm = []
    mask = []
    for i in range(atms.shape[1]):
        sel = np.array(sels[:,i],bool)
        n, m = normalizeAtm(atms[:,i],sel,stable)
        norm.append(n)
        mask.append(m)
        # sel = np.array(sels[:,i],bool)*stable
        # if (sel.sum() > 40) and (sels[:,i].sum() > 100): 
        #     m =  np.mean(atms[:,i][sel])
        #     atms[:,i] /= m
        #     norm.append(m)
        #     mask.append(True)
        # else:
        #     norm.append(1)
        #     mask.append(False)
    return np.array(norm),numpy.array(mask)


# DEPRECATED
def normalizeAtm(atms,sel,stable,min_stable=40,min_sel=100):
    sel_st = sel*stable
    if (sel_st.sum() > min_stable) and (sel.sum() > min_sel):
        m =  np.mean(atms[sel_st])
        atms /= m
        mask = True
    else:
        m = 1
        mask = False
    return m, mask

# DEPRECATED
def normalizeAtm_stable(atms,ff,resp,sel,stable):
    # Warning: Raw atms should be used here!
    sel_st = sel*stable
    if (sel_st.sum() > 20) and (sel.sum() > 100):
        m =  np.mean((resp*ff)[sel_st]) / np.mean((resp*atms)[sel_st])
        atms /= m
        mask = True
    else:
        m = 1
        mask = False
    return m, mask


def findMode(x,window=1.2e-5,init=[0,7e-5],maxiter=20,minerr=1e-9):
    sel = (x>init[0])*(x<init[1])
    ss = window/2
    m0,s=medsig(x[sel])
    for i in range(maxiter):
        sel = (x>m0-ss)*(x<m0+ss)
        m,s=medsig(x[sel])
        #print m0,m,s
        if abs(m-m0) < minerr:
            return m,s
        else:
            m0 = m
    print("WARNING: Maximum iterations reaached finding distribution mode")
    return m,s
    
def medsig(x):
    xsi = numpy.array(x).argsort()
    n = len(x)
    m = x[xsi[int(n/2)]]
    q25 = x[xsi[int(n/4)]]
    q75 = x[xsi[int(3.*n/4.)]]
    s = 0.741*(q75-q25)
    return m,s
