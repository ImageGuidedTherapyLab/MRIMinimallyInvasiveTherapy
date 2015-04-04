import sqlite3
import time
import dicom
import os
#conn = sqlite3.connect('/data/fuentes/DICOM/ctkDICOM.sql')
conn = sqlite3.connect('./ctkDICOM.sql')
# build second database
os.system('rm ./pydicomtags.sql')
tagsconn = sqlite3.connect('./pydicomtags.sql')

# build view for tags
# verify tag info
# select * from dicom_header h1 where h1.dicomkey = "(0018, 0081)";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 0080)";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 0082)";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 0087)";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 0095)";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 0088)";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 0050)";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 0083)";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 0086)";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 1314)";
# select * from dicom_header h1 where h1.dicomkey = "(0008, 0020)";
# select * from dicom_header h1 where h1.dicomkey = "(0008, 103e)";
# select * from dicom_header h1 where h1.name like "%slab%";
# select * from dicom_header h1 where h1.dicomkey = "(0018, 9104)";
# select * from dicom_header h1 where h1.name like "%pulse%";

# tagsconn.execute("""
# create view adniinfo as 
#    select h00.name,h00.value,
#           h01.name,h01.value,
#           h02.name,h02.value,
#           h03.name,h03.value,
#           h04.name,h04.value,
#           h05.name,h05.value,
#           h06.name,h06.value,
#           h07.name,h07.value,
#           h08.name,h08.value,
#           h09.name,h09.value,
#           h10.name,h10.value,
#           h11.name,h11.value,
#           h12.name,h12.value,
#           h13.name,h13.value,
#           h14.name,h14.value,
#           h15.name,h15.value,
#           h16.name,h16.value,
#           h17.name,h17.value,
#           h18.name,h18.value
#    from dicom_header h00
#    join dicom_header h01 on h00.SeriesInstanceUID=h01.SeriesInstanceUID
#    join dicom_header h02 on h00.SeriesInstanceUID=h02.SeriesInstanceUID
#    join dicom_header h03 on h00.SeriesInstanceUID=h03.SeriesInstanceUID
#    join dicom_header h04 on h00.SeriesInstanceUID=h04.SeriesInstanceUID
#    join dicom_header h05 on h00.SeriesInstanceUID=h05.SeriesInstanceUID
#    join dicom_header h06 on h00.SeriesInstanceUID=h06.SeriesInstanceUID
#    join dicom_header h07 on h00.SeriesInstanceUID=h07.SeriesInstanceUID
#    join dicom_header h08 on h00.SeriesInstanceUID=h08.SeriesInstanceUID
#    join dicom_header h09 on h00.SeriesInstanceUID=h09.SeriesInstanceUID
#    join dicom_header h10 on h00.SeriesInstanceUID=h10.SeriesInstanceUID
#    join dicom_header h11 on h00.SeriesInstanceUID=h11.SeriesInstanceUID
#    join dicom_header h12 on h00.SeriesInstanceUID=h12.SeriesInstanceUID
#    join dicom_header h13 on h00.SeriesInstanceUID=h13.SeriesInstanceUID
#    join dicom_header h14 on h00.SeriesInstanceUID=h14.SeriesInstanceUID
#    join dicom_header h15 on h00.SeriesInstanceUID=h15.SeriesInstanceUID
#    join dicom_header h16 on h00.SeriesInstanceUID=h16.SeriesInstanceUID
#    join dicom_header h17 on h00.SeriesInstanceUID=h17.SeriesInstanceUID
#    join dicom_header h18 on h00.SeriesInstanceUID=h18.SeriesInstanceUID
#    where h00.dicomkey = "(0008, 103e)" 
#    and   h01.dicomkey = "(0018, 0086)" 
#    and   h02.dicomkey = "(0008, 0020)"
#    and   h03.dicomkey = "(0018, 0081)"
#    and   h04.dicomkey = "(0018, 0080)"
#    and   h05.dicomkey = "(0018, 0082)"
#    and   h06.dicomkey = "(0018, 0087)"
#    and   h07.dicomkey = "(0018, 0095)"
#    and   h08.dicomkey = "(0018, 0088)"
#    and   h09.dicomkey = "(0018, 0050)"
#    and   h10.dicomkey = "(0018, 0083)"
#    and   h11.dicomkey = "(0018, 1314)"
#    and   h12.dicomkey = "(0019, 109c)"
#    and   h13.dicomkey = "(0018, 0094)"
#    and   h14.dicomkey = "(0018, 1310)"
#    and   h15.dicomkey = "(0028, 0011)"
#    and   h16.dicomkey = "(0028, 0010)"
#    and   h17.dicomkey = "(0021, 1056)"
#    and   h18.dicomkey = "(0021, 1057)"
#    ;
# """)

## # drop table dicom_header; 
## # delete from dicom_header; 
## #conn.execute('create table if not exists dicom_header (SOPInstanceUID varchar(64) not null, dicomkey varchar(32) not null ,name varchar(64), value varchar(512), primary key (SOPInstanceUID,dicomkey) )')
## conn.execute('create table if not exists dicom_header (SeriesInstanceUID VARCHAR(64) NOT NULL , dicomkey varchar(32) not null ,name varchar(64), value varchar(512), PRIMARY KEY (SeriesInstanceUID,dicomkey) )' )
## 
## #fileIDList = [ (seriesUID,filename) for (seriesUID,filename) in conn.execute('select SeriesInstanceUID,Filename from Images limit 5')]
## 
## # only update files that are "DaysBackOld"
## DaysBack = 2;
## # convert all to seconds
## # 60 sec * 60 min * 24hr = number of secs in day
## datecutoff = time.strftime("%Y-%m-%d",time.localtime(time.time() - 60.*60.*24.*DaysBack )) 
## # select by modality and date
## fileIDList = [ (seriesUID,filename) for (seriesUID,filename) in conn.execute('select im.SeriesInstanceUID,im.Filename from Images im join Series se on im.SeriesInstanceUID=se.SeriesInstanceUID where DateTime(im.InsertTimeStamp) > "%s" and se.Modality like "%%mr%%" ' % datecutoff )]
## 
## 
# select files 
## fileIDList = [ (filename,seriesUID) for (filename,seriesUID) in conn.execute('select im.Filename,se.SeriesInstanceUID from ( (Studies sd join Series se on se.StudyInstanceUID=sd.StudyInstanceUID) join Images im on se.SeriesInstanceUID=im.SeriesInstanceUID ) where sd.StudyDate = "2013-07-17" group by se.SeriesDescription ;')]

tagsconn.execute('create table if not exists TagCache (SOPInstanceUID VARCHAR NOT NULL , Tag varchar not null ,Name varchar, Value varchar, PRIMARY KEY (SOPInstanceUID,Tag) )' )
fileIDList = [ (filename,sopUID) for (filename,sopUID) in conn.execute('select Filename,SOPInstanceUID from Images ')]

# build database
errlogfileHandle = file('err.txt'  ,'w')
for (filename,sopUID) in fileIDList:
  print 'Processing %s' % filename
  dcm=dicom.read_file(filename);
  # loop over all keys
  for dcmkey in dcm.keys() :
      try: 
        # catch key exceptions
        name=dcm[dcmkey].name;
        value=dcm[dcmkey].value;
        tableentry=(unicode(str(sopUID)),unicode(str(dcmkey)),unicode(str(name)),unicode(str(value)))
        # insert and ignore duplicate entries
      ## except sqlite3.IntegrityError as inst:
      ## except UnicodeDecodeError as inst:
      ## except ValueError as inst:
      # catch key exceptions
      except Exception as inst:
        name='NameException'
        value='ValueException'
        tableentry=(unicode(str(sopUID)),unicode(str(dcmkey)),unicode(str(name)),unicode(str(value)))
        errlogfileHandle.write("%s," % inst )
        errlogfileHandle.write('%s,%s\n' %  ( unicode(str(sopUID)),unicode(str(dcmkey)) ) )
      finally:
        # this is always executed
        tagsconn.execute('insert or ignore into TagCache (SOPInstanceUID,Tag,Name,Value) values (?,?,?,?);' , tableentry)
        errlogfileHandle.flush()
  tagsconn.commit()

errlogfileHandle.close()







## USEFUL COMMANDS
## ---------------
## select count(*) from dicom_header where Name="Slice Thickness" and cast(value as number)<2;
## select value from dicom_header where name="Slice Thickness" group by value order by value;
## select im.SeriesInstanceUID,se.SeriesDescription,h.value from dicom_header h join Images im on im.SOPInstanceUID=h.SOPInstanceUID join Series se on se.SeriesInstanceUID=im.SeriesInstanceUID where h.Name="Slice Thickness" and cast(h.value as number)<=2 group by im.SeriesInstanceUID;

## DROP view IF EXISTS 'tmpsubset' ;
## create view tmpsubset as select * from dicom_header where value = 012345 and name = "Patient ID" ; 
## select sb.name,sb.value,hd.name,hd.value,im.SeriesInstanceUID,im.Filename from dicom_header hd join tmpsubset sb on sb.SOPInstanceUID=hd.SOPInstanceUID join Images im on hd.SOPInstanceUID=im.SOPInstanceUID where hd.Name like "%slice%" and cast(hd.value as number) <=3 ; 

## drop view WS_PACS;
## create view WS_PACS as select st.StudyInstanceUID,se.SeriesInstanceUID,im.FileName from images im join Series se on se.SeriesInstanceUID=im.SeriesInstanceUID join Studies st on st.StudyInstanceUID=se.StudyInstanceUID;

## sqlite> select sr.SeriesDescription,hd.value from series sr join Images im on sr.SeriesInstanceUID = im.SeriesInstanceUID join dicom_header hd on im.SOPInstanceUID = hd.SOPInstanceUID where sr.SeriesDescription like '%ax%' and  hd.dicomkey = '(0020, 0037)' limit 4;
## sqlite> select sr.SeriesDescription,hd.value from series sr join Images im on sr.SeriesInstanceUID = im.SeriesInstanceUID join dicom_header hd on im.SOPInstanceUID = hd.SOPInstanceUID where sr.SeriesDescription like '%sag%' and  hd.dicomkey = '(0020, 0037)' limit 4;
## sqlite> select sr.SeriesDescription,hd.value from series sr join Images im on sr.SeriesInstanceUID = im.SeriesInstanceUID join dicom_header hd on im.SOPInstanceUID = hd.SOPInstanceUID where sr.SeriesDescription like '%cor%' and  hd.dicomkey = '(0020, 0037)' limit 4;

# select hd.value,tr.name,tr.value,te.name,te.value from dicom_header hd join echotime te on te.SeriesInstanceUID=hd.SeriesInstanceUID join reptime tr on hd.SeriesInstanceUID=tr.SeriesInstanceUID where hd.value like '%T1%' or hd.value like '%T2%' ;
