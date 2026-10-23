"""Download the public processed export; refuse empty or changed responses."""
import argparse
import hashlib
import urllib.request
from pathlib import Path

URL='https://api.figshare.com/v2/file/download/34464122'
SHA='e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830'
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args()
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==SHA:
        print('Existing checksum-verified dataset retained.')
    else:
        temporary=path.with_suffix('.partial')
        urllib.request.urlretrieve(URL,temporary)
        if hashlib.sha256(temporary.read_bytes()).hexdigest()!=SHA:
            temporary.unlink(missing_ok=True)
            raise ValueError('Download checksum mismatch; original destination was not changed')
        temporary.replace(path)
        print('Downloaded and verified processed dataset.')
