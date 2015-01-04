struct GetRequest {
    int64_t block_id;
    int offset;
    int length;
    char* data;
    int timeout_ms;
};

struct PutRequest {
    int64_t block_id;
    std::map<int, int64_t> server_addrs;
    const char* data;
    int length;
    char* code;
    int timeout_ms;
};

struct Part {
    int length;
    int offset;
    char* data;
};

struct PartsRequest {
    int64_t block_id;
    std::map<int, int64_t> server_addrs;
    std::map<int, int64_t> avail_server_addrs;
    std::map<int, Part> todo_parts;
    std::map<uint32_t, Part> finish_parts;
    int timeout_ms;
};

struct PartRequest {

};

int RSClient::Get(GetRequest* request) {
    //1 check request parameter
    //2 split the data buffer into several parts, where a part belongs to one replica
    //3 query meta handle to get block replicas location
    //4 query part handle to get parts data
}

int PartHandle::GetParts(PartsRequest* request) {
    //for each todo parts do GetPart();
}

int RSClient::Put(PutRequest* request) {
    //1 check request parameter
    //2 encode request on code handle
    //2 make PartsRequest from request
    //3 part_handle_->PutParts(request);
    //4 time and conditional wait
}

int PartHandle::PutParts(PartsRequest* request) {
    //for each todo parts to make PartRequest and TryPutPart()
}

int PartHandle::TryPutPart() {
    //control max running put part request
    //do PutPart() : it is implemented in the data server
}

int PartHandle::OnFinishPutPart() {
    if all done, FinishRequest
    else TryPutPart();
}

