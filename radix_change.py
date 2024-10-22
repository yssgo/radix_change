import re, fractions, sys

# from addsep.py
import itertools

_sepchars={
    2:"_",
    8:"_",
    10:",",
    16:"_",
}
_sepchars2={
    2:"_",
    8:"_",
    10:"",
    16:"_",
}

_partsizes={
    2:4,
    8:4,
    10:3,
    16:2
}


def partition(s,size):
    def _split_seq(s,size):
        it=iter(list(s))
        item=list(itertools.islice(it,size))
        while item:
            yield item
            item=list(itertools.islice(it,size))
    return [''.join(_) for _ in _split_seq(s,size)]


def addsep2(s,sep,partsize, fromend=True, padchar=""):
    if not fromend:
        ret=list(partition(s,partsize))
    else:
        ret = list(partition(s[::-1],partsize))
    if padchar != "":
        last = ret[-1]
        rest = partsize - len(last)
        if rest >=1:
            ret[-1] = last + padchar*rest
    if not fromend:
        ret=sep.join(ret)   
    else:
        ret = sep.join(ret)[::-1]
    return ret  


def addsep(s,sep1=",",sep2="",partsize=3, padchar=""):
    if "." in s:
        ipart,fpart= s.split(".")
        ipart2= addsep2(ipart,sep1,partsize, fromend=True, padchar=padchar)
        fpart2= addsep2(fpart,sep2,partsize, fromend=False, padchar=padchar)
        return ipart2+"."+fpart2
    else:
        ipart = s
        ipart2 = addsep2(ipart,sep1, partsize, fromend=True, padchar=padchar)
        return ipart2


def addseparator(s,radix=10,decsep2="", partsize=0, padchar=''):
    if radix in _sepchars:
        sep=_sepchars[radix]
    else:
        sep="_"
    if radix in _sepchars2:
        if radix !=10:
            sep2=_sepchars2[radix]
        else:
            sep2=decsep2
    else:
        sep2="_"
    if partsize == 0:
        if radix in _partsizes:
            psize=_partsizes[radix]
        else:
            psize=3
    else:
        psize = partsize
    return addsep(s,sep,sep2,psize,padchar=padchar)
#end : from addsep.py

global _MAXBINFRAC
global _rptstart, _rptstop
global _destip, _DESTFP

def isBinInfinte(sn10):     
    ipart, fpart = sn10.split(".")
    fr=fractions.Fraction(int(fpart), 10**(len(fpart)))
    num = fr.denominator
    
    binstr = ""
    
    while True:
        q = num // 2
        r = num % 2
        binstr =str(r)+binstr
        if q==0:
            break
        num=q
    
    return binstr.count('1') != 1
        
def getRepeatedInfo():
    return _DESTFP, _rptstart, _rptstop

_MAXBINFRAC = 23
def setBinPrecision(precision):
    global _MAXBINFRAC
    _MAXBINFRAC = precision
_digitchars="0123456789"+\
"ABCDEFGHIJKLMNOPQRSTUVWWXYZ"
_bintable= {
    16: [f"{_:04b}" for _ in range(16)],
    8: [f"{_:03b}" for _ in range(8)]
}

class RadixChange:
    @classmethod
    def _makedestfp_not10(self,srcfp, from_radix, to_radix):
        global _DESTFP
        if from_radix==to_radix:
            _DESTFP = srcfp
            return _DESTFP
        _DESTFP=""
        if from_radix==2:
            partlens ={8:3, 16:4}
            partlen = partlens[to_radix]
            while len(srcfp) % partlen != 0:
                srcfp += "0"
            for i in range(0, len(srcfp),partlen):
              part = srcfp[i:i+partlen]
              partval =0
              for i, e in enumerate(part[::-1]):
                  partval += _digitchars.find(e) * (2 ** i)
              _DESTFP+=_digitchars[partval]
        else:
            for i in range(len(srcfp)):
                digit_val = _digitchars.find(srcfp[i])
                _DESTFP += _bintable[from_radix][digit_val]
            if to_radix==2:
                return _DESTFP
            partlens ={8:3, 16:4}
            partlen = partlens[to_radix]
            while len(_DESTFP)%partlen!=0:
                _DESTFP+="0"
            srcfp=_DESTFP[:]
            _DESTFP=""
            for i in range(0,len(srcfp),partlen):
                part=srcfp[i:i+partlen]
                _DESTFP+=_bintable[to_radix][int(part,2)]
        return _DESTFP
    @classmethod
    def _makedestfp_to10(self, srcfp,from_radix):
        binf = 0
        for i, e in enumerate(srcfp):
            binf += _digitchars.find(e) * (from_radix**-(i+1))
        return str(binf).partition('.')[2]
    @classmethod
    def _makedestfp_from10(self, srcfp, from_radix, to_radix, show=True):
        global _rptstart,_rptstop,_DESTFP,CALCEDLEN
        _rptstart=-1
        _rptstop=-1
        _DESTFP=""
        srcfp_val = int(srcfp,from_radix)
        srcfp_val_len = len(str(srcfp_val))
        num = srcfp_val
        rpt=[]
        
        pos=-1
        while True:
            pos +=1
            if _rptstart==-1 and num in rpt:
                _rptstart=rpt.index(num)
                _rptstop=pos
            else:
                rpt.append(num)
            n1 = num
            n2 = n1 * to_radix
            n1s=f"{n1:0{srcfp_val_len+1}d}"
            n1ip=n1s[0]
            n1fp=n1s[1:].rstrip('0')
            if n1fp=="":
                n1fp="0"
            n2s=f"{n2:0{srcfp_val_len+1}d}"
            n2ip=n2s[0]
            n2fp=n2s[1:].rstrip("0")
            if n2fp=="":
                n2fp="0"
            if(show): print(f"{n1ip}.{n1fp} * {to_radix} = {n2ip}.{n2fp}", end="")
            #if show and to_radix!=10:
            #   print(f" ( {n2ip}={_digitchars[int(n2ip)]}({to_radix}) )", end="")
            if(show): print()
            q = n2 // (10**(srcfp_val_len))
            r = n2 % (10**(srcfp_val_len))
            _DESTFP += _digitchars[q]
            if r == 0 or len(_DESTFP) == _MAXBINFRAC:
                CALCEDLEN=len(_DESTFP)
                break
            num = r
        return _DESTFP
    def __init__():
        pass
    @classmethod
    def radix_change(self, sn, from_radix, to_radix, decsep2="_", show=True, nosep=False,partsize=0, padchar=''):
        global _rptstart, _rptstop
        global destip, _DESTFP,CALCEDLEN
        CALCEDLEN= 0
        _rptstart=-1; _rptstop=-1
        #sn = str(sn)
    
        if '.' in sn:
            srcip,srcfp = sn.strip().upper().split('.')
        else:
            srcip, srcfp = sn.strip().upper(), ""
        if from_radix == to_radix:
            destip = srcip
            destfp = srcfp
            if destfp!="":
                outs=destip+"."+destfp
            else:
                outs=destip
            
            if not nosep: outs=addseparator(outs, to_radix,decsep2=decsep2, partsize=partsize, padchar=padchar)
            return outs

        srcip_val = int(srcip, from_radix)

        destip = ""
        if srcip_val == 0:
            destip="0"
        else:
            n = srcip_val
            while True:
                q = n // to_radix
                r = n % to_radix
                if(show): print(f"{n} / {to_radix} = {q} ... {r}", end="")
                #if show and to_radix!=10:
                #   print(f" ( {r}={_digitchars[r]}({to_radix}) )", end="")
                if(show): print()
                destip = _digitchars[r] + destip
                if q == 0:
                    break
                n = q
        if srcfp !="":
            if from_radix == 10:
                destfp = self._makedestfp_from10(srcfp, from_radix, to_radix)
            elif to_radix == 10:
                destfp = self._makedestfp_to10(srcfp, from_radix)
            else:
                destfp = self._makedestfp_not10(srcfp, from_radix, to_radix)
        else:
            destfp=""
            
        if destfp!="":
            outs=destip+"."+destfp
        else:
            outs=destip
        if not nosep:            
            outs=addseparator(outs, to_radix, decsep2=decsep2, partsize=partsize, padchar=padchar)
        return outs

def stripsep(sn, seps=r"[_',]"):
    return re.sub(seps,"",sn)


def tellrep(sn, from_radix):
    if from_radix!=10:
        return ""
    if "." in sn:
        if isBinInfinte(sn):
            print("무한 소수임")
        else:
            print("유한 소수임")
    else:
        return ""
    if _rptstart!=-1:        
        def rplen2ko(rpstr):
            rplen = len(rpstr)
            if 1<=rplen<=10:
                return ["한","두","세","네","다섯","여섯","일곱","여덟","아홉","열"]
            else:
                return str(rplen)
        rpstr = _DESTFP[_rptstart:_rptstop]
        rpstrko = rplen2ko(rpstr)
        def numjosa(num):
            if str(num)[-1].upper() in '013678LMN':
                return "이"
            else:
                return "가"
        _DESTFP[_rptstart:_rptstop]
        print("소수점 아래",
         (_rptstart+1), "번째 자리에서부터", _DESTFP[_rptstart:_rptstop], f"{numjosa(_DESTFP[_rptstop-1])} 무한 반복됨({rpstrko} 자리)")
    else:
        if CALCEDLEN >= 1:
            print(f"소수점 아래 {CALCEDLEN} 자리까지 반복되는 수 없음.")
        else:
            print("소수점 아래 반복되는 수 없음.")

if __name__=="__main__":
    from_radix = -1
    while from_radix not in [2,8,10,16]:
        try:
            choice = input("입력 진수(2/8/10/16). Enter=취소: ")
            if choice == '':
                sys.exit(0)
            from_radix = int(choice)
        except ValueError:
            print("잘못된 입력 진수입니다.")
            continue;
        if from_radix not in [2,8,10,16]:
            print("잘못된 입력 진수입니다.")
    to_radix = -1
    while to_radix not in set([2,8,10,16])-set([from_radix]):
        to_radix = -1
        try:
            choice = input(f"출력 진수({'/'.join(map(str,sorted(set([2,8,10,16])-set([from_radix]))))}). Enter=취소: ")
            if choice == '':
                sys.exit(0)
            to_radix = int(choice)
        except ValueError: 
            print("ValueError")
            print("출력 진수에 대한 입력이 잘못되었습니다.")
            print("f{to_radix=}, {type(to_radix)=}")
            continue;
        if to_radix not in set([2,8,10,16])-set([from_radix]):
            print("잘못된 출력 진수입니다.")
            print("f{to_radix=}, {type(to_radix)=}")
    
    passed = 0
    while not passed:
        passed = 0
        sn = input(str(from_radix)+"진수(양수만 가능)를 입력하십시오(Enter=취소): ")
        if sn=='':
            sys.exit(0)
        sn = re.sub("[_,' \t]","",sn)
        if all(x in '.'+"0123456789ABCDEFGHIJLMNOPQRSTUVWXYZ"[:from_radix] for x in sn):
            if sn.count('.') in [0,1]:
                passed = 1
    
    while True:
        mbf=input(f"소수점 이하 최대 몇 자리까지 구할까요? (기본값: {_MAXBINFRAC}) Enter=취소: ")
        if mbf=="":
            break
        try:
            _MAXBINFRAC = int(mbf)
            break
        except ValueError:
            print("잘못된 입력입니다.")
            continue
    print()
    
    
    
    if to_radix==2:
        binsn = RadixChange.radix_change(sn,from_radix,2,partsize=4, padchar='0')
    if to_radix==8:
        octsn = RadixChange.radix_change(sn,from_radix,8,padchar='0')
    if to_radix==10:
        decsn = RadixChange.radix_change(sn,from_radix,10) 
    if to_radix==16:
        hexsn = RadixChange.radix_change(sn,from_radix,16,padchar='0')        
        
    if to_radix == 2:
        print(binsn,"(2)")        
        print(addseparator(stripsep(binsn), 2, partsize=3, padchar='0'),'(2)')
        tellrep(sn, from_radix)
        print()
    if to_radix==8:
        print(octsn,"(8)")
        tellrep(sn, from_radix)
        print()
    if to_radix==10:
        print(decsn,"(10)")
        print()

    if to_radix==16:
        print(hexsn,"(16)")
        tellrep(sn, from_radix)
        print()
    
    # if to_radix ==2:
        # print(addseparator(stripsep(binsn), radix=2, partsize=4, padchar='0'),'(2)')
    # if to_radix == 8:
        # binsn = RadixChange.radix_change(sn,from_radix,2,partsize=3, padchar='0', show=0)
        # print(addseparator(stripsep(binsn), radix=2, partsize=3, padchar='0'),'(2)')        
    # if to_radix == 8:
        # print(octsn,'(8)')
    # if to_radix == 10:
        # print(decsn,'(10)')
    # if to_radix == 16:
        # print(hexsn,'(16)')
    
    if to_radix == 2:
        print(stripsep(binsn),"(2)")    
    if to_radix == 8:
        print(stripsep(octsn),"(8)")
    if to_radix == 10:
        print(stripsep(decsn),"(10)")
    # if to_radix==10:
        # print(stripsep(decsn,r"[_]"),"(10)")
    if to_radix == 16:
        print(stripsep(hexsn),"(16)")    
        