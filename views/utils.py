#! /sw/arch/bin/python
#
# *************************************************************
#
# $Source: $
# $Revision: $                                                                 
# $State: $                                                                     
# $Date: $                                                      
# $Author: $  
#
# $Log: $
#
#
# *************************************************************

def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result



def pheno4j():
    neo=get_neo4j()
    result = neo.run("MATCH (a:Person) return a.personId as personId ")
    s='\n'.join([ "%s" % (record["personId"]) for record in result ])
    return s



def individuals_update(external_ids):
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    users_db=get_db(app.config['DB_NAME_USERS'])
    def f(eid):
        p=patients_db.patients.find_one({'external_id':eid},{'_id':False})
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        if '_id' in p: del p['_id']
        return p
    new_individuals=[f(eid) for eid in external_ids]
    old_individuals=users_db.users.find_one({'user':session['user']}).get('individuals',[])
    old_individuals=[ind for ind in old_individuals if ind['external_id'] not in external_ids]
    individuals=new_individuals+old_individuals
    users_db.users.update_one({'user':session['user']},{'$set':{'individuals':individuals}})
    return individuals


def get_feature_venn(patient):
    s="""
    MATCH (p:Person)-[:PersonToObservedTerm]->(t:Term)--(g:Gene)
    WHERE p.personId='%s'
    RETURN t.termId, t.name, g.gene_id, g.gene_name
    """ % patient
    print(s)
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)
    data = []
    for r in result:
        data.append({
            'hpo_id': r['t.termId'],
            'hpo_term': r['t.name'],
            'gene_id': r['g.gene_id'],
            'gene_name': r['g.gene_name']
        })
    hpo_terms=[(k,v,) for k, v, in dict([(x['hpo_id'],x['hpo_term'],) for x in data]).items()]
    hpo_gene=dict()
    for x in data:
        hpo_gene[x['hpo_id']]=hpo_gene.get(x['hpo_id'],[])+[x['gene_name']]
    genes = {}
    feature_combo = []
    feature_venn = []
    print "get combinatorics of features to draw venn diagram"
    for i in range(len(hpo_terms[:5])):
        feature_combo.extend(itertools.combinations(range(len(hpo_terms)), i+1))
    print 'calculate Venn diagram'
    for combo in feature_combo:
        # construct features_venn key
        #venn_ind += 1
        dic_key = [hpo_terms[i][1] for i in combo]
        for ind in range(len(combo)):
            if ind == 0:
                x=hpo_terms[combo[ind]][0]
                feature_venn.append({'key': dic_key, 'value':list(frozenset(hpo_gene.get(x,"")))})
            else:
                tem = feature_venn[-1]['value']
                feature_venn[-1]['value'] = list(frozenset(feature_venn[-1]['value']) & frozenset(hpo_gene[hpo_terms[combo[ind]][0]]))
    return feature_venn

def patient_variants():
    # add known gene and retnet gene labels, and re-calculate pubmed_score
    for mm in ['rare_variants','homozygous_variants','compound_het_variants']:
        for v in patient.__dict__[mm]:
            if 'canonical_gene_name_upper' not in v: v['canonical_gene_name_upper']=v['Gene']
            gene=v['canonical_gene_name_upper']
            pubmed_key = '_'.join([gene,patient.get('pubmed_key','')])
            gene_info[gene]=dict()
            if gene in known_genes: 
                gene_info[gene]['known']=True
                pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if gene not in RETNET: continue
            gene_info[gene]['disease'] = RETNET[gene]['disease']
            gene_info[gene]['omim'] = RETNET[gene]['omim']
            gene_info[gene]['mode'] = RETNET[gene]['mode']
            pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if mm != 'rare_variants' or ('d' in gene_info[gene]['mode'] and mm == 'rare_variants') :
                pubmedbatch[pubmed_key] = max(100,pubmedbatch[pubmed_key])
                if gene=='DRAM2':
                    print pubmed_key
                    print pubmedbatch[pubmed_key]
            if 'het_samples' not in v: print(v)
            for s in v['het_samples']:
                if v['HET_COUNT'] < 10:
                    individuals[s]=individuals.get(s,[])+[v]


def exomiser(individual):
    patient_hpo_terms=lookups.get_patient_hpo(hpo_db, patient_db, individual, ancestors=False)
    patient_hpo_terms = dict([(hpo['id'][0],{'id':hpo['id'][0],'name':hpo['name'][0], 'is_a':hpo.get('is_a',[])}) for hpo in patient_hpo_terms])
    patient_hpo_ids=patient_hpo_terms.keys()
    x['exomiser']=[]
    for g in list(set(x['genes'])):
        r=db.ensembl_entrez.find_one({'Ensembl Gene ID':g})
        if not r or not r['EntrezGene ID']: continue
        x['entrezgeneid']=r['EntrezGene ID']
        #url='http://localhost:8085/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        url='http://monarch-exomiser-prod.monarchinitiative.org/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        print(url)
        r=requests.get(url)
        if isinstance(r.json(),list):
            x['exomiser']+=r.json()[0]['results']
        else:
            x['exomiser']+=r.json()['results']
    if len(x['exomiser'])<1: x['exomiser']=[{'score':-1}]
    exomiser_scores=[xx['score'] for xx in x['exomiser']]
    i=exomiser_scores.index(max(exomiser_scores))
    x['exomiser']=x['exomiser'][i]



def patient_variants():
    # add known gene and retnet gene labels, and re-calculate pubmed_score
    for mm in ['rare_variants','homozygous_variants','compound_het_variants']:
        for v in patient.__dict__[mm]:
            if 'canonical_gene_name_upper' not in v: v['canonical_gene_name_upper']=v['Gene']
            gene=v['canonical_gene_name_upper']
            pubmed_key = '_'.join([gene,patient.get('pubmed_key','')])
            gene_info[gene]=dict()
            if gene in known_genes: 
                gene_info[gene]['known']=True
                pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if gene not in RETNET: continue
            gene_info[gene]['disease'] = RETNET[gene]['disease']
            gene_info[gene]['omim'] = RETNET[gene]['omim']
            gene_info[gene]['mode'] = RETNET[gene]['mode']
            pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if mm != 'rare_variants' or ('d' in gene_info[gene]['mode'] and mm == 'rare_variants') :
                pubmedbatch[pubmed_key] = max(100,pubmedbatch[pubmed_key])
                if gene=='DRAM2':
                    print pubmed_key
                    print pubmedbatch[pubmed_key]
            if 'het_samples' not in v: print(v)
            for s in v['het_samples']:
                if v['HET_COUNT'] < 10:
                    individuals[s]=individuals.get(s,[])+[v]


def get_hpo_gene(hpo_ids):
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    hpo_terms = [(i, hpo_db.hpo.find_one({'id':i})['name'][0]) for i in hpo_ids]
    # this has missing HPO ids. see IRDC_batch2_OXF_3001 and #HP:0000593
    hpo_gene=dict()
    for hpo_id,hpo_term, in hpo_terms:
        hpo_gene[hpo_id] = []
        for gene_name in [x['Gene-Name'] for x in hpo_db.ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.find({'HPO-ID':hpo_id},{'Gene-Name':1,'_id':0})]:
            #gene_hpo[gene_name]=gene_hpo.get(gene_name,[])+[{'hpo_id':hpo_id,'hpo_term':hpo_term}]
            hpo_gene[hpo_id]=hpo_gene.get(hpo_id,[])+[gene_name]
    for k in hpo_gene: hpo_gene[k]=list(frozenset(list(hpo_gene[k])))
    return hpo_gene


def exomiser(individual):
    patient_hpo_terms=lookups.get_patient_hpo(hpo_db, patient_db, individual, ancestors=False)
    patient_hpo_terms = dict([(hpo['id'][0],{'id':hpo['id'][0],'name':hpo['name'][0], 'is_a':hpo.get('is_a',[])}) for hpo in patient_hpo_terms])
    patient_hpo_ids=patient_hpo_terms.keys()
    x['exomiser']=[]
    for g in list(set(x['genes'])):
        r=db.ensembl_entrez.find_one({'Ensembl Gene ID':g})
        if not r or not r['EntrezGene ID']: continue
        x['entrezgeneid']=r['EntrezGene ID']
        #url='http://localhost:8085/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        url='http://monarch-exomiser-prod.monarchinitiative.org/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        print(url)
        r=requests.get(url)
        if isinstance(r.json(),list):
            x['exomiser']+=r.json()[0]['results']
        else:
            x['exomiser']+=r.json()['results']
    if len(x['exomiser'])<1: x['exomiser']=[{'score':-1}]
    exomiser_scores=[xx['score'] for xx in x['exomiser']]
    i=exomiser_scores.index(max(exomiser_scores))
    x['exomiser']=x['exomiser'][i]




def individuals_update(external_ids):
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    users_db=get_db(app.config['DB_NAME_USERS'])
    def f(eid):
        p=patients_db.patients.find_one({'external_id':eid},{'_id':False})
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        if '_id' in p: del p['_id']
        return p
    new_individuals=[f(eid) for eid in external_ids]
    old_individuals=users_db.users.find_one({'user':session['user']}).get('individuals',[])
    old_individuals=[ind for ind in old_individuals if ind['external_id'] not in external_ids]
    individuals=new_individuals+old_individuals
    users_db.users.update_one({'user':session['user']},{'$set':{'individuals':individuals}})
    return individuals


def get_feature_venn(patient):
    s="""
    MATCH (p:Person)-[:PersonToObservedTerm]->(t:Term)--(g:Gene)
    WHERE p.personId='%s'
    RETURN t.termId, t.name, g.gene_id, g.gene_name
    """ % patient
    print(s)
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)
    data = []
    for r in result:
        data.append({
            'hpo_id': r['t.termId'],
            'hpo_term': r['t.name'],
            'gene_id': r['g.gene_id'],
            'gene_name': r['g.gene_name']
        })
    hpo_terms=[(k,v,) for k, v, in dict([(x['hpo_id'],x['hpo_term'],) for x in data]).items()]
    hpo_gene=dict()
    for x in data:
        hpo_gene[x['hpo_id']]=hpo_gene.get(x['hpo_id'],[])+[x['gene_name']]
    genes = {}
    feature_combo = []
    feature_venn = []
    print "get combinatorics of features to draw venn diagram"
    for i in range(len(hpo_terms[:5])):
        feature_combo.extend(itertools.combinations(range(len(hpo_terms)), i+1))
    print 'calculate Venn diagram'
    for combo in feature_combo:
        # construct features_venn key
        #venn_ind += 1
        dic_key = [hpo_terms[i][1] for i in combo]
        for ind in range(len(combo)):
            if ind == 0:
                x=hpo_terms[combo[ind]][0]
                feature_venn.append({'key': dic_key, 'value':list(frozenset(hpo_gene.get(x,"")))})
            else:
                tem = feature_venn[-1]['value']
                feature_venn[-1]['value'] = list(frozenset(feature_venn[-1]['value']) & frozenset(hpo_gene[hpo_terms[combo[ind]][0]]))
    return feature_venn


def load_patient(individual,auth,pubmed_key,hpo='HP:0000001'):
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    db = get_db()
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    patient_id=individual
    patient={u'features': {u'observed': u'yes', u'type': u'phenotype', u'id': hpo}, 'clinicalStatus': {u'clinicalStatus': u'affected'}, u'ethnicity': {u'maternal_ethnicity': [], u'paternal_ethnicity': []}, u'family_history': {}, u'disorders': [], u'life_status': u'alive', u'reporter': u'', u'genes': [], u'prenatal_perinatal_phenotype': {u'prenatal_phenotype': [], u'negative_prenatal_phenotype': []}, u'prenatal_perinatal_history': {u'twinNumber': u''}, u'sex': u'U', u'solved': {u'status': u'unsolved'}}
    eid=patient_id
    if p: patient.update(p)
    #patient_hpo_terms=','.join([f['id'] for f in patient['features'] if f['observed']=='yes'])
    gene_counter=Counter([var['canonical_gene_name_upper'] for var in patient.rare_variants])
    for var in patient['rare_variants']: var['gene_count']=gene_counter[var['canonical_gene_name_upper']]
    patient["pubmedbatch_status"]=0
    pubmed_key="blindness-macula-macular-pigmentosa-retina-retinal-retinitis-stargardt"
    patient["pubmed_key"]=pubmed_key
    #db.patients.update({'external_id':patient_id}, patient, upsert=True)



def pheno4j():
    neo=get_neo4j()
    result = neo.run("MATCH (a:Person) return a.personId as personId ")
    s='\n'.join([ "%s" % (record["personId"]) for record in result ])
    return s


def individuals_update(external_ids):
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    users_db=get_db(app.config['DB_NAME_USERS'])
    def f(eid):
        p=patients_db.patients.find_one({'external_id':eid},{'_id':False})
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        if '_id' in p: del p['_id']
        return p
    new_individuals=[f(eid) for eid in external_ids]
    old_individuals=users_db.users.find_one({'user':session['user']}).get('individuals',[])
    old_individuals=[ind for ind in old_individuals if ind['external_id'] not in external_ids]
    individuals=new_individuals+old_individuals
    users_db.users.update_one({'user':session['user']},{'$set':{'individuals':individuals}})
    return individuals



def patient_variants():
    # add known gene and retnet gene labels, and re-calculate pubmed_score
    for mm in ['rare_variants','homozygous_variants','compound_het_variants']:
        for v in patient.__dict__[mm]:
            if 'canonical_gene_name_upper' not in v: v['canonical_gene_name_upper']=v['Gene']
            gene=v['canonical_gene_name_upper']
            pubmed_key = '_'.join([gene,patient.get('pubmed_key','')])
            gene_info[gene]=dict()
            if gene in known_genes: 
                gene_info[gene]['known']=True
                pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if gene not in RETNET: continue
            gene_info[gene]['disease'] = RETNET[gene]['disease']
            gene_info[gene]['omim'] = RETNET[gene]['omim']
            gene_info[gene]['mode'] = RETNET[gene]['mode']
            pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if mm != 'rare_variants' or ('d' in gene_info[gene]['mode'] and mm == 'rare_variants') :
                pubmedbatch[pubmed_key] = max(100,pubmedbatch[pubmed_key])
                if gene=='DRAM2':
                    print pubmed_key
                    print pubmedbatch[pubmed_key]
            if 'het_samples' not in v: print(v)
            for s in v['het_samples']:
                if v['HET_COUNT'] < 10:
                    individuals[s]=individuals.get(s,[])+[v]


def exomiser(individual):
    patient_hpo_terms=lookups.get_patient_hpo(hpo_db, patient_db, individual, ancestors=False)
    patient_hpo_terms = dict([(hpo['id'][0],{'id':hpo['id'][0],'name':hpo['name'][0], 'is_a':hpo.get('is_a',[])}) for hpo in patient_hpo_terms])
    patient_hpo_ids=patient_hpo_terms.keys()
    x['exomiser']=[]
    for g in list(set(x['genes'])):
        r=db.ensembl_entrez.find_one({'Ensembl Gene ID':g})
        if not r or not r['EntrezGene ID']: continue
        x['entrezgeneid']=r['EntrezGene ID']
        #url='http://localhost:8085/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        url='http://monarch-exomiser-prod.monarchinitiative.org/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        print(url)
        r=requests.get(url)
        if isinstance(r.json(),list):
            x['exomiser']+=r.json()[0]['results']
        else:
            x['exomiser']+=r.json()['results']
    if len(x['exomiser'])<1: x['exomiser']=[{'score':-1}]
    exomiser_scores=[xx['score'] for xx in x['exomiser']]
    i=exomiser_scores.index(max(exomiser_scores))
    x['exomiser']=x['exomiser'][i]

def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def setup_neo4j_driver(host, port, password):
    default_password = 'neo4j' # Travis will use a fresh Neo4j, with the default password.
    local_password = password
    uri = "bolt://"+host+":"+str(port)
    # Try local_password.
    try:
        driver = GraphDatabase.driver(uri, auth=basic_auth("neo4j", local_password))
        return driver
    except:
        pass
    # Try default_password.
    # Password handling from https://github.com/robinedwards/django-neomodel
    driver = GraphDatabase.driver(uri, auth=basic_auth("neo4j", default_password))
    with driver.session() as neo4j_session:
        try:
            result = neo4j_session.run("MATCH (a:Person) WHERE a.name = {name} RETURN a", {"name": "Crick"})

        except CypherError as ce:
            if 'The credentials you provided were valid, but must be changed before you can use this instance' in str(ce):
                neo4j_session.run("CALL dbms.changePassword({password})", {'password': local_password})
                print("New database with no password set, setting password to '", local_password, "'.")
                neo4j_session.close()
            else:
                raise ce
    return driver

def create_demo_user(neo4j_session):
    results = neo4j_session.run("MATCH (u:User {user : 'demo'}) RETURN u")
    result = results.single()
    if not result:
        print("Adding user 'demo' to the neo4j database.")
        neo4j_session.run("CREATE (a:User {user: {username}, argon_password: {hash}})", {"username": "demo", "hash": argon2.hash("demo123")}) 



def get_feature_venn(patient):
    s="""
    MATCH (p:Person)-[:PersonToObservedTerm]->(t:Term)--(g:Gene)
    WHERE p.personId='%s'
    RETURN t.termId, t.name, g.gene_id, g.gene_name
    """ % patient
    print(s)
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)
    data = []
    for r in result:
        data.append({
            'hpo_id': r['t.termId'],
            'hpo_term': r['t.name'],
            'gene_id': r['g.gene_id'],
            'gene_name': r['g.gene_name']
        })
    hpo_terms=[(k,v,) for k, v, in dict([(x['hpo_id'],x['hpo_term'],) for x in data]).items()]
    hpo_gene=dict()
    for x in data:
        hpo_gene[x['hpo_id']]=hpo_gene.get(x['hpo_id'],[])+[x['gene_name']]
    genes = {}
    feature_combo = []
    feature_venn = []
    print "get combinatorics of features to draw venn diagram"
    for i in range(len(hpo_terms[:5])):
        feature_combo.extend(itertools.combinations(range(len(hpo_terms)), i+1))
    print 'calculate Venn diagram'
    for combo in feature_combo:
        # construct features_venn key
        #venn_ind += 1
        dic_key = [hpo_terms[i][1] for i in combo]
        for ind in range(len(combo)):
            if ind == 0:
                x=hpo_terms[combo[ind]][0]
                feature_venn.append({'key': dic_key, 'value':list(frozenset(hpo_gene.get(x,"")))})
            else:
                tem = feature_venn[-1]['value']
                feature_venn[-1]['value'] = list(frozenset(feature_venn[-1]['value']) & frozenset(hpo_gene[hpo_terms[combo[ind]][0]]))
    return feature_venn

def individuals_update(external_ids):
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    users_db=get_db(app.config['DB_NAME_USERS'])
    def f(eid):
        p=patients_db.patients.find_one({'external_id':eid},{'_id':False})
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        if '_id' in p: del p['_id']
        return p
    new_individuals=[f(eid) for eid in external_ids]
    old_individuals=users_db.users.find_one({'user':session['user']}).get('individuals',[])
    old_individuals=[ind for ind in old_individuals if ind['external_id'] not in external_ids]
    individuals=new_individuals+old_individuals
    users_db.users.update_one({'user':session['user']},{'$set':{'individuals':individuals}})
    return individuals


def load_patient(individual,auth,pubmed_key,hpo='HP:0000001'):
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    db = get_db()
    patient_db=get_db(app.config['DB_NAME_PATIENTS'])
    patient_id=individual
    patient={u'features': {u'observed': u'yes', u'type': u'phenotype', u'id': hpo}, 'clinicalStatus': {u'clinicalStatus': u'affected'}, u'ethnicity': {u'maternal_ethnicity': [], u'paternal_ethnicity': []}, u'family_history': {}, u'disorders': [], u'life_status': u'alive', u'reporter': u'', u'genes': [], u'prenatal_perinatal_phenotype': {u'prenatal_phenotype': [], u'negative_prenatal_phenotype': []}, u'prenatal_perinatal_history': {u'twinNumber': u''}, u'sex': u'U', u'solved': {u'status': u'unsolved'}}
    eid=patient_id
    if p: patient.update(p)
    #patient_hpo_terms=','.join([f['id'] for f in patient['features'] if f['observed']=='yes'])
    gene_counter=Counter([var['canonical_gene_name_upper'] for var in patient.rare_variants])
    for var in patient['rare_variants']: var['gene_count']=gene_counter[var['canonical_gene_name_upper']]
    patient["pubmedbatch_status"]=0
    pubmed_key="blindness-macula-macular-pigmentosa-retina-retinal-retinitis-stargardt"
    patient["pubmed_key"]=pubmed_key
    #db.patients.update({'external_id':patient_id}, patient, upsert=True)


def individuals_update(external_ids):
    patients_db=get_db(app.config['DB_NAME_PATIENTS'])
    users_db=get_db(app.config['DB_NAME_USERS'])
    def f(eid):
        p=patients_db.patients.find_one({'external_id':eid},{'_id':False})
        print p['external_id']
        p['features']=[f for f in p.get('features',[]) if f['observed']=='yes']
        if 'solved' in p:
            if 'gene' in p['solved']:
                p['solved']=[p['solved']['gene']]
            else:
                p['solved']=[]
        else: p['solved']=[]
        if 'genes' in p: p['genes']=[x['gene'] for x in p['genes'] if 'gene' in x]
        else: p['genes']=[]
        p['genes']=list(frozenset(p['genes']+p['solved']))
        p2=get_db().patients.find_one({'external_id':p['external_id']},{'rare_homozygous_variants_count':1,'rare_compound_hets_count':1, 'rare_variants_count':1,'total_variant_count':1})
        if not p2: return p
        p['rare_homozygous_variants_count']=p2.get('rare_homozygous_variants_count','')
        p['rare_compound_hets_count']=p2.get('rare_compound_hets_count','')
        p['rare_variants_count']=p2.get('rare_variants_count','')
        p['total_variant_count']=p2.get('total_variant_count','')
        #p['all_variants_count']=get_db().patients.find_one({'external_id':p['external_id']},{'_id':0,'all_variants_count':1})['all_variants_count']
        #db.cache.find_one({"key" : "%s_blindness,macula,macular,retina,retinal,retinitis,stargardt_" % })
        if '_id' in p: del p['_id']
        return p
    new_individuals=[f(eid) for eid in external_ids]
    old_individuals=users_db.users.find_one({'user':session['user']}).get('individuals',[])
    old_individuals=[ind for ind in old_individuals if ind['external_id'] not in external_ids]
    individuals=new_individuals+old_individuals
    users_db.users.update_one({'user':session['user']},{'$set':{'individuals':individuals}})
    return individuals


def phenogenon(hpo_id,lit_genes,omim_genes,recessive_genes,dominant_genes,cache=True):
    cache_db=get_db('cache')
    temp=cache_db.phenogenon_cache.find_one({'hpo_id':hpo_id})
    if temp and cache:
        lit_genes.extend(temp['lit_genes'])
        omim_genes.extend(temp['omim_genes'])
        recessive_genes.extend(temp['recessive_genes'])
        dominant_genes.extend(temp['dominant_genes'])
        return
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    db=get_db()
    def f(r):
        g=db.genes.find_one({'gene_name_upper':r['Gene-Name'].upper()},{'_id':0})
        if not g: return
        phenogenon=db.gene_hpo.find_one({'gene_id':g['gene_id']})
        if not phenogenon: return
        het=phenogenon.get('het',{}).get(hpo_id,{})
        hom_comp=phenogenon.get('hom_comp',{}).get(hpo_id,{})
        if 'data' in het: del het['data']
        if 'data' in hom_comp: del hom_comp['data']
        g['phenogenon']={ 'het':het, 'hom_comp': hom_comp}
        return g
    lit_genes=[f(r) for r in hpo_db.hpo_gene.find({'HPO-ID':hpo_id})]
    lit_genes=[lg for lg in lit_genes if lg]
    omim_genes.extend(map(lambda x: x['gene_id'], lit_genes))
    phenogenon=db.hpo_gene.find_one({'hpo_id':hpo_id})
    if phenogenon: phenogenon=phenogenon['data']['unrelated']
    else: phenogenon={'recessive':[],'dominant':[]}
    recessive_genes.extend([{'gene_id':x['gene_id'],'gene_name':db.genes.find_one({'gene_id':x['gene_id']})['gene_name'],'p_val':x['p_val'],'known':x['gene_id'] in omim_genes} for x in phenogenon['recessive']])
    dominant_genes.extend([{'gene_id':x['gene_id'],'gene_name':db.genes.find_one({'gene_id':x['gene_id']})['gene_name'],'p_val':x['p_val'], 'known':x['gene_id'] in omim_genes} for x in phenogenon['dominant']])
    #print({'hpo_id':hpo_id,'dominant_genes':dominant_genes,'recessive_genes':recessive_genes,'omim_genes':omim_genes,'lit_genes':lit_genes})
    cache_db.phenogenon_cache.insert_one({'hpo_id':hpo_id,'dominant_genes':dominant_genes,'recessive_genes':recessive_genes,'omim_genes':omim_genes,'lit_genes':lit_genes})

def skat(hpo_id):
    db=get_db()
    skat_genes=db.skat.find({'HPO':hpo_id},{'_id':False})
    skat_genes=[g for g in skat_genes if g['FisherPvalue']<0.05 and g['SKATO']<0.005]
    for g in skat_genes:
        pli=get_db('exac').pli.find_one({'gene':g['Symbol']})
        if str(g['OddsRatio'])=='inf': g['OddsRatio']=9999999
        if pli:
            g['pli']=pli['pLI'] 
        else:
            g['pli']=-1
    return skat_genes


def get_feature_venn(patient):
    s="""
    MATCH (p:Person)-[:PersonToObservedTerm]->(t:Term)--(g:Gene)
    WHERE p.personId='%s'
    RETURN t.termId, t.name, g.gene_id, g.gene_name
    """ % patient
    print(s)
    with neo4j_driver.session() as neo4j_session:
        result=neo4j_session.run(s)
    data = []
    for r in result:
        data.append({
            'hpo_id': r['t.termId'],
            'hpo_term': r['t.name'],
            'gene_id': r['g.gene_id'],
            'gene_name': r['g.gene_name']
        })
    hpo_terms=[(k,v,) for k, v, in dict([(x['hpo_id'],x['hpo_term'],) for x in data]).items()]
    hpo_gene=dict()
    for x in data:
        hpo_gene[x['hpo_id']]=hpo_gene.get(x['hpo_id'],[])+[x['gene_name']]
    genes = {}
    feature_combo = []
    feature_venn = []
    print "get combinatorics of features to draw venn diagram"
    for i in range(len(hpo_terms[:5])):
        feature_combo.extend(itertools.combinations(range(len(hpo_terms)), i+1))
    print 'calculate Venn diagram'
    for combo in feature_combo:
        # construct features_venn key
        #venn_ind += 1
        dic_key = [hpo_terms[i][1] for i in combo]
        for ind in range(len(combo)):
            if ind == 0:
                x=hpo_terms[combo[ind]][0]
                feature_venn.append({'key': dic_key, 'value':list(frozenset(hpo_gene.get(x,"")))})
            else:
                tem = feature_venn[-1]['value']
                feature_venn[-1]['value'] = list(frozenset(feature_venn[-1]['value']) & frozenset(hpo_gene[hpo_terms[combo[ind]][0]]))
    return feature_venn



def patient_variants():
    # add known gene and retnet gene labels, and re-calculate pubmed_score
    for mm in ['rare_variants','homozygous_variants','compound_het_variants']:
        for v in patient.__dict__[mm]:
            if 'canonical_gene_name_upper' not in v: v['canonical_gene_name_upper']=v['Gene']
            gene=v['canonical_gene_name_upper']
            pubmed_key = '_'.join([gene,patient.get('pubmed_key','')])
            gene_info[gene]=dict()
            if gene in known_genes: 
                gene_info[gene]['known']=True
                pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if gene not in RETNET: continue
            gene_info[gene]['disease'] = RETNET[gene]['disease']
            gene_info[gene]['omim'] = RETNET[gene]['omim']
            gene_info[gene]['mode'] = RETNET[gene]['mode']
            pubmedbatch[pubmed_key] = max(1,pubmedbatch.get('pubmed_key',0))
            if mm != 'rare_variants' or ('d' in gene_info[gene]['mode'] and mm == 'rare_variants') :
                pubmedbatch[pubmed_key] = max(100,pubmedbatch[pubmed_key])
                if gene=='DRAM2':
                    print pubmed_key
                    print pubmedbatch[pubmed_key]
            if 'het_samples' not in v: print(v)
            for s in v['het_samples']:
                if v['HET_COUNT'] < 10:
                    individuals[s]=individuals.get(s,[])+[v]


def get_hpo_gene(hpo_ids):
    hpo_db=get_db(app.config['DB_NAME_HPO'])
    hpo_terms = [(i, hpo_db.hpo.find_one({'id':i})['name'][0]) for i in hpo_ids]
    # this has missing HPO ids. see IRDC_batch2_OXF_3001 and #HP:0000593
    hpo_gene=dict()
    for hpo_id,hpo_term, in hpo_terms:
        hpo_gene[hpo_id] = []
        for gene_name in [x['Gene-Name'] for x in hpo_db.ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.find({'HPO-ID':hpo_id},{'Gene-Name':1,'_id':0})]:
            #gene_hpo[gene_name]=gene_hpo.get(gene_name,[])+[{'hpo_id':hpo_id,'hpo_term':hpo_term}]
            hpo_gene[hpo_id]=hpo_gene.get(hpo_id,[])+[gene_name]
    for k in hpo_gene: hpo_gene[k]=list(frozenset(list(hpo_gene[k])))
    return hpo_gene


def exomiser(individual):
    patient_hpo_terms=lookups.get_patient_hpo(hpo_db, patient_db, individual, ancestors=False)
    patient_hpo_terms = dict([(hpo['id'][0],{'id':hpo['id'][0],'name':hpo['name'][0], 'is_a':hpo.get('is_a',[])}) for hpo in patient_hpo_terms])
    patient_hpo_ids=patient_hpo_terms.keys()
    x['exomiser']=[]
    for g in list(set(x['genes'])):
        r=db.ensembl_entrez.find_one({'Ensembl Gene ID':g})
        if not r or not r['EntrezGene ID']: continue
        x['entrezgeneid']=r['EntrezGene ID']
        #url='http://localhost:8085/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        url='http://monarch-exomiser-prod.monarchinitiative.org/exomiser/api/prioritise/?phenotypes=%s&prioritiser=hiphive&genes=%s&prioritiser-params=human,mouse,fish'%(','.join(patient_hpo_terms.keys()), x['entrezgeneid'])
        print(url)
        r=requests.get(url)
        if isinstance(r.json(),list):
            x['exomiser']+=r.json()[0]['results']
        else:
            x['exomiser']+=r.json()['results']
    if len(x['exomiser'])<1: x['exomiser']=[{'score':-1}]
    exomiser_scores=[xx['score'] for xx in x['exomiser']]
    i=exomiser_scores.index(max(exomiser_scores))
    x['exomiser']=x['exomiser'][i]


'''
defs
'''
def hide_id_for_demo(data):
    if not data: return
    for k,v in data['patients'].items():
        # hide hpo
        v['hpo'] = ['hidden']
        # hide variants
        v['variants'] = ['hidden_'+hashlib.sha224(i).hexdigest()[:6] for i in v['variants']]
        # hide p_id
        new_p = 'hidden_'+hashlib.sha224(k).hexdigest()[:6]
        data['patients'][new_p] = data['patients'].pop(k)
    for k1,v1 in data['data'].items():
        for k2,v2 in v1['p'].items():
            v1['p'][k2] = ['hidden_'+hashlib.sha224(i).hexdigest()[:6] for i in v2]
    for k,v in data['variants'].items():
        new_v = 'hidden_'+hashlib.sha224(k).hexdigest()[:6]
        data['variants'][new_v] = data['variants'].pop(k)




# For use in easy_install.sh to set up a demo user.
if __name__ == '__main__':
    uri = sys.argv[1] 
    password = sys.argv[2] 
    driver = GraphDatabase.driver(uri, auth=basic_auth("neo4j", password))
    with driver.session() as neo4j_session:  
        create_demo_user(neo4j_session)
