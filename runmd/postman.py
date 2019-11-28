#/usr/bin/python

import os
import sys
import smtplib
# For guessing MIME type based on file name extension
import mimetypes


from email import encoders
from email.mime.audio import MIMEAudio
from email.mime.base import MIMEBase
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import getpass


def send_mail(subject, text, to=None, cc=None, files=None):
    """
    :type cc: list of str
    :type to: list of str
    :type files: list of str
    """

    assert isinstance(to, list) or to is None
    assert isinstance(cc, list) or cc is None
    assert isinstance(files, list) or files is None
    send_from = "bionmr.server@gmail.com"



    usersmail = {
      "root" : "sergei.a.izmailov@gmail.com",
      "admin" : "sergei.a.izmailov@gmail.com",
      "sergei":"sergei.a.izmailov@gmail.com",
      "seva" : "sevastyan@rabdano.ru",
      "ivan" : "podkorytov@rambler.ru",
      "egorov" : "egorov@nmr.phys.spbu.ru",
      "kerstin" : "kerstin.kaempf@gmx.de",
      "ubovich" : "ubovich.milosh@yandex.ru",
      "olga" : "acerlaetum@ya.ru"
      }

    username=getpass.getuser()

    if to is None: to=[]
    if cc is None: cc=[]

    to = [ usersmail[username] ] + to
    cc = cc + [ usersmail["admin"] ]


    msg = MIMEMultipart()

    msg["From"] = send_from
    msg["To"] = ", ".join(to)
    msg["Cc"] = ", ".join(cc)
    msg["Subject"] = subject
    msg['Reply-to'] = usersmail["admin"]

    msg.attach(MIMEText(text))

    for path in files or []:
        ctype, encoding = mimetypes.guess_type(path)
        if ctype is None or encoding is not None:
            # No guess could be made, or the file is encoded (compressed), so
            # use a generic bag-of-bits type.
            ctype = 'application/octet-stream'
        maintype, subtype = ctype.split('/', 1)

        if not os.path.isfile(path):
            print >> sys.stderr, "WARNING: File `%s` is not regular file or not exists, so not attached. " % path
            continue

        if maintype == 'text':
            fp = open(path)
            # Note: we should handle calculating the charset
            inner = MIMEText(fp.read(), _subtype=subtype)
            fp.close()
        elif maintype == 'image':
            fp = open(path, 'rb')
            inner = MIMEImage(fp.read(), _subtype=subtype)
            fp.close()
        elif maintype == 'audio':
            fp = open(path, 'rb')
            inner = MIMEAudio(fp.read(), _subtype=subtype)
            fp.close()
        else:
            fp = open(path, 'rb')
            inner = MIMEBase(maintype, subtype)
            inner.set_payload(fp.read())
            fp.close()
            # Encode the payload using Base64
            encoders.encode_base64(inner)
        # Set the filename parameter
        inner.add_header('Content-Disposition', 'attachment', filename=path)
        msg.attach(inner)

    smtp = smtplib.SMTP('smtp.gmail.com')
    smtp.ehlo()
    smtp.starttls()
    smtp.login(send_from, "snparknmqsqyzpcb")
    smtp.sendmail(send_from, to+cc, msg.as_string())
    smtp.close()
