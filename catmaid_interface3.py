# catmaid_interface: A library for interacting with CATMAID projects through python
# Started by Tom Kazimiers (1/2013)
# Adapted by Albert Cardona (1/2013)
# Casey Schneider-Mizell (8/2015)

import urllib.request, urllib.parse, urllib.error
import ssl
import base64
import http.cookiejar
import json

# Class for a basic CATMAID login/fetch/POST/GET with authentication
class Connection:
    def __init__(self, server, authname, authpassword, username, password):
        self.server = server
        self.authname = authname
        self.authpassword = authpassword
        self.username = username
        self.password = password
        self.cookies = http.cookiejar.CookieJar()
        self.opener = urllib.request.build_opener(urllib.request.HTTPRedirectHandler(), urllib.request.HTTPCookieProcessor(self.cookies))
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        self.ctx = ctx

    def djangourl(self, path):
        """ Expects the path to lead with a slash '/'. """
        return self.server + path

    def auth(self, request):
        if self.authname:
            loginstr = self.authname + ':' + self.authpassword
            loginstr.replace('\n', '')
            base64string = base64.b64encode(bytes('%s:%s' % (self.authname, self.authpassword),'utf-8') )
            request.add_header(b'Authorization', b'Basic ' + base64string)

    def login(self):
        url = self.djangourl("/accounts/login")
        opts = {
            b'name': self.username,
            b'pwd': self.password
        }
        data = urllib.parse.urlencode(opts)
        request = urllib.request.Request(url, bytes(data,'utf-8') )
        self.auth(request)
        response = urllib.request.urlopen(request, context=self.ctx)
        self.cookies.extract_cookies(response, request)
        return response.read()

    def fetch(self, url, post=None):
        """ Requires the url to connect to and the variables for POST, if any, in a dictionary. """
        if post:
            request = urllib.request.Request(url, bytes(post,'utf-8'))
        else:
            request = urllib.request.Request(url)

        self.auth(request)
        return self.opener.open(request).read()
        
    def fetchGET(self, url, data ):
        full_url = url + '?' + data
        response = self.fetch(full_url)
        if not response:
            return
        r = json.loads(response.decode('utf-8'))
        if type(r) == dict and 'error' in r:
            print(("ERROR:", r['error']))
        else:
            return r        

    def fetchJSON(self, url, post=None):
        response = self.fetch(url, post=post)
        if not response:
            return
        r = json.loads(response.decode('utf-8'))
        if type(r) == dict and 'error' in r:
            print(("ERROR:", r['error']))
        else:
            return r

# Given a list of skeleton_ids and a connection, gets the JSON file for all of them
def write_skeletons_from_list( c, id_list, project_id ):
    for id in id_list:
        sk = skeleton_object(c, project_id, id)
        if not not sk:
            f_nodes = open( 'sk_' + str(id) + '.json' ,'w')
            sk.append( id )
            sk.append( get_neuron_name(c, project_id, id) )
            presyn_weight = post_synaptic_count(c,project_id, [conn[1] for conn in sk[1] if conn[2] == 0])
            if presyn_weight!=None:
                sk.append([ [conn[0], len(conn[1]['postsynaptic_to'])] for conn in presyn_weight])
            else:
                sk.append([])
            json.dump(sk,f_nodes)
            f_nodes.close()    

# Add one annotation to a list of skeleton IDs.
def add_annotation( c, project_id, annotation_str, id_list ):
    anno_url = c.djangourl('/' + str(project_id) + "/annotations/add")
    anno_post = dict()
    anno_post["annotations[0]"] = annotation_str
    for i, id in enumerate(id_list):
        anno_post["skeleton_ids[" + str(i) + "]"] = id
    anno_post_url = urllib.parse.urlencode( anno_post )
    d = c.fetch( anno_url, post=anno_post_url )
    return d
    
# get_neuron_name: Given a skeleton ID, fetch the neuron name
def get_neuron_name( c, project_id, skeleton_id ):
    url = c.djangourl('/%s/skeleton/%s/neuronname' % (project_id, skeleton_id) )
    d = c.fetchJSON(url)
    return d['neuronname']
                
# Fetch JSON for one skeleton
def skeleton_object( c, project_id, skeleton_id):
    url = c.djangourl('/%s/%s/1/1/compact-skeleton' % (project_id, skeleton_id))
    d = c.fetchJSON(url)
    return d

# Count how many postsynaptic targets are associated with a list of connectors               
def post_synaptic_count( c, project_id, connector_list ):
    url = c.djangourl('/%s/connector/skeletons' % (project_id))
    opts = {}
    for ind, id in enumerate(connector_list):
        opts[ 'connector_ids['+str(ind)+']' ] = id
    data = urllib.parse.urlencode(opts)
    returned_info = c.fetchJSON(url,data)
    return returned_info

# Parse a file to a list of skeleton ids.        
def import_ids(filename):
    f = open(filename,'r')
    id_list = []
    for line in f:
        if int(line) == -1:
            f.close()
            return []
        else:
            id_list.append(int(line))
    f.close()
    return id_list

# Grow a list of skeleton ids along the connectivity network     
def increase_id_list( c, id_list_base, min_pre, min_post, hops, project_id):
    if len(id_list_base)==0:
        return []
    
    opts = {'n_circles': hops,
        'min_pre': min_pre,
        'min_post': min_post}
    for ind, id in enumerate(id_list_base):
        opts[ 'skeleton_ids['+str(ind)+']' ] = id

    url = c.djangourl('/%s/graph/circlesofhell' % (project_id) )
    data = urllib.parse.urlencode(opts)
    
    returned_info = c.fetchJSON(url,data)
    id_list = returned_info[0]
    
    for id in id_list_base:
        id_list.append( id )
    id_list = list(set(id_list))
    return id_list

# Given an annotation id, pull all skeletons with that annotation
def get_ids_from_annotation(c, project_id, annotation_id_list):
    url = c.djangourl('/%s/neuron/query-by-annotations' % (project_id) )
    id_list = []
    for anno_id in annotation_id_list:
        opts = {'neuron_query_by_annotation' : anno_id }
        data = urllib.parse.urlencode(opts)
        returned_info = c.fetchJSON(url,data)
        ids_in_returned_info = [ entry['skeleton_ids'][0] for entry in returned_info['entities'] ]
        id_list = id_list + ids_in_returned_info
    id_list = list(set(id_list))
    return id_list

# Gets the dict of annotation_id to annotation string for a project                
def make_annotation_dict(c, project_id ):
    anno_url = c.djangourl('/%s/annotations/list' % (project_id))
    all_annotations = c.fetchJSON( anno_url )
    anno_dict = { item['name']:item['id'] for item in all_annotations['annotations'] }
    return anno_dict
    
# Fetch ids of all skeletons with more nodes than min_nodes
def fetch_nodes_by_count( c, project_id, min_nodes=1000 ):
    url = c.djangourl('/%s/skeleton/list' % (project_id) )
    opts = {'nodecount_gt': min_nodes}
    data = urllib.parse.urlencode(opts)
    return_list = c.fetchGET( url, data )
    return return_list