classdef Pushbullet < handle
    %Pushbullet Connect to pushbullet.com and send notes/links/files to the
    %connected smartphone as a push notification
    %
    %   Usage:
    %   p = Pushbullet('abcdefghijk123456') (apikey)
    %   p.pushNote([],'Matlab Pushbullet Test','This is the message') --
    %   sends a push note to all connected devices
    %   p.pushNote('abhgzt12123','Matlab Pushbullet Test','This is the
    %   message') --send only to the specific device with the device_iden
    %   p.load_devices()  --show all devices + device_idens
    %   p.pushFile(device_iden, file_name, file_type, file_url) to push a
    %   file which has already been uploaded
    %
    %   Copyright 2016, Jens Brauer, https://github.com/jensb89
    
    properties
        HOST = 'https://api.pushbullet.com/v2'
        DEVICES_URL = '/devices'
        CONTACTS_URL = '/contacts'
        ME_URL = '/users/me'
        PUSH_URL = '/pushes'
        UPLOAD_REQUEST_URL = '/upload-request'
        SUBSCRIPTIONS_URL = '/subscriptions'
        ApiKey
        Devices
        Subscriptions
    end
    
    methods
        
        function self = Pushbullet(apikey)
            self.ApiKey = apikey;
            if verLessThan('MATLAB','8.5')
                warning(['You are using a Matlab version prior to 2015a',...
                        'Matlab-Pushbullet might not work with this version.',...
                        'Update Matlab or use an older version of Matlab-Pushbullet',...
                        'from here: https://github.com/jensb89/Matlab-Pushbullet/releases/']);
            end
        end
        
        function load_devices(self)
            % Get a list of devices
            options =  weboptions('KeyName','Access-Token','KeyValue',self.ApiKey);
            output = webread([self.HOST,self.DEVICES_URL],options);
            self.Devices = output.devices;
            for i=1:length(self.Devices)
                if self.Devices{i}.active && isfield(self.Devices{i},'nickname')
                sprintf('%s : %s', self.Devices{i}.nickname, self.Devices{i}.iden)
                end
            end
        end
        
        function load_subscriptions(self)
            % Get a list of devices
            options =  weboptions('KeyName','Access-Token','KeyValue',self.ApiKey);
            output = webread([self.HOST,self.SUBSCRIPTIONS_URL],options);
            self.Subscriptions = output.subscriptions;
            for i=1:length(self.Subscriptions)
                if self.Subscriptions{i}.active && isfield(self.Subscriptions{i}.channel,'name')
                sprintf('%s : %s : %s', self.Subscriptions{i}.channel.name, self.Subscriptions{i}.channel.iden, self.Subscriptions{i}.channel.tag)
                end
            end
        end      
        
        
        function output = pushNote(self, device_iden, title, message)
            % Push a note
            % https://docs.pushbullet.com/v2/pushes
            % Arguments:
            % device_iden -- iden of device to push to
            % title -- a title for the note
            % body -- the body of the note

            data = struct('body',message,...
                          'device_iden',device_iden,...
                          'title',title,...
                          'type','note');
            
            if isempty(device_iden)
                data = rmfield(data,'device_iden');
            end
            
            output = push(self, data);
        end
        
     function output = pushNote_to_Channel(self, channel_tag, title, message)
            % Push a note
            % https://docs.pushbullet.com/v2/pushes
            % Arguments:
            % channel_tag -- tag of channel to push to
            % title -- a title for the note
            % body -- the body of the note

            data = struct('body',message,...
                          'channel_tag',channel_tag,...
                          'title',title,...
                          'type','note');
            
            if isempty(channel_tag)
                data = rmfield(data,'channel_tag');
            end
            
            output = push(self, data);
        end
 
                      
        function output = pushFile(self, device_iden, file_name, file_type, file_url)
            % Push a picture
            % https://docs.pushbullet.com/v2/pushes
            % Arguments:
            % device_iden -- iden of device to push to
            % file_name -- the name for the file
            % file_type -- The MIME type of the file (for example
            % 'image/png')
            % file_url -- the url of the file

            data = struct('type','file',...
                    'device_iden',device_iden,...
                    'file_name',file_name,...
                    'file_type',file_type,...
                    'file_url',file_url);
                
            if isempty(device_iden)
                data = rmfield(data,'device_iden'); %delete device_iden in data -> push to all connected devices
            end

            output = push(self, data);
        end
        
        function output = pushLink(self, device_iden, title, message, url)
            % Push a link
            % https://docs.pushbullet.com/v2/pushes
            % Arguments:
            % device_iden -- iden of device to push to
            % title -- a title for the note
            % body -- the body of the note

            data = struct('type', 'note',...
                'device_iden',device_iden,...
                'title', title,...
                'body', message,...
                'url', url);
            
            if isempty(device_iden)
                data = rmfield(data,'device_iden'); %delete device_iden in data -> push to all connected devices
            end
            
            output = push(self, data);
        end
            
        function output = push(self, data)
            % Perform the POST Request
            options = weboptions('KeyName','Access-Token',...
                                'KeyValue',self.ApiKey,...
                                'MediaType','application/json');
            output = webwrite([self.HOST,self.PUSH_URL],data,options);
        end
        
        % HELPER FUNCTIONS
        function device_iden = get_device_iden_from_nickname(self, nickname)
            if isempty(self.Devices)
                load_devices(self)
            end
            if ~isempty(self.Devices)
                for i=1:length(self.Devices)
                    if strcmp(self.Devices{i}.nickname, nickname)
                        device_iden = self.Devices{i}.iden;
                    end
                end
            end
        end       
    
    end
    
end

