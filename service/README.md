# systemctl service to autostart webpage

* paths are hard coded
* If SELinux runs this may fail:
  * check `journalctl`


* enable service (see below for failures):
```bash
ln -s /opt/GCMAMBOX/service/gcmambox.service /etc/systemd/system/

systemctl daemon-reload
systemctl enable gcmambox.service
systemctl start gcmambox.service
systemctl status gcmambox.service
```
